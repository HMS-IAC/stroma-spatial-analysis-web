/**
 * Authors: Antoine A. Ruzette, Simon F. Nørrelykke
 * Date: 2024-03-20
 *
 * This script annotates stromal regions based on a threshold of the Fibronectin marker, then 
 * calculates the signed distance between cells and their closest stromal border.
 * It also supports scanning through different threshold combinations and calculates the
 * Pearson's correlation coefficient between the signed distance and the intensity of the marker.
 * 
 * Released under the MIT License (see LICENSE file)
 */

import qupath.lib.objects.PathObject
import static qupath.lib.gui.scripting.QPEx.*
import java.nio.file.Files
import java.nio.file.Paths
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation

// Iterate over threshold combinations
def sigmas = (10..11).step(1).collect()
def fn_thresholds = (4700..5000).step(300).collect()


//def ker_thresholds = (0..2000000).step(100000).collect()
//def marker_thresholds = (0..500000).step(25000).collect()

println('Starting stroma annotation...')

sigmas.each { sigma ->
    fn_thresholds.each { fn_threshold ->

        selectObjects { p -> p.getPathClass() == getPathClass("Stroma") && p.isAnnotation() }
        def selectedObjects = QP.getSelectedObjects()

        if (!selectedObjects.isEmpty()) {
            removeObjects(selectedObjects, false)
            println("Selected Stroma annotations deleted.")
        } else {
            println("No Stroma annotations selected for deletion.")
        }

        resetSelection()
        selectAnnotations()

        println("Parameters: Sigma: ${sigma}; FN_568_threshold: ${fn_threshold}")
        
        def pixelClassifierString = """
        {
          "pixel_classifier_type": "OpenCVPixelClassifier",
          "metadata": {
            "inputPadding": 0,
            "inputResolution": {
              "pixelWidth": {
                "value": 1.2924389465034254,
                "unit": "µm"
              },
              "pixelHeight": {
                "value": 1.292467583464835,
                "unit": "µm"
              },
              "zSpacing": {
                "value": 1.0,
                "unit": "z-slice"
              },
              "timeUnit": "SECONDS",
              "timepoints": []
            },
            "inputWidth": 512,
            "inputHeight": 512,
            "inputNumChannels": 3,
            "outputType": "CLASSIFICATION",
            "outputChannels": [],
            "classificationLabels": {
              "0": {},
              "1": {
                "name": "Stroma",
                "color": [
                  150,
                  200,
                  150
                ]
              }
            }
          },
          "op": {
            "type": "data.op.channels",
            "colorTransforms": [
              {
                "channelName": "TRITC FN"
              }
            ],
            "op": {
              "type": "op.core.sequential",
              "ops": [
                {
                  "type": "op.filters.gaussian",
                  "sigmaX": ${sigma},
                  "sigmaY": ${sigma}
                },
                {
                  "type": "op.threshold.constant",
                  "thresholds": [
                    ${fn_threshold}
                  ]
                }
              ]
            }
          }
        }
        """

        def pixelClassifierFilePath = "/classifiers/stroma_annotation/stroma_annotator_sigma=${sigma}_fn-thresh=${fn_threshold}.json"
        new File(pixelClassifierFilePath).text = pixelClassifierString

        createAnnotationsFromPixelClassifier(pixelClassifierFilePath, 0.0, 0.0)
        println('Done!')

        println('Starting calculation of signed distance...')
        detectionToAnnotationDistancesSigned(false)
        println('Done!')

        println("Starting calculation of Pearson's correlation coefficients...")

        ArrayList<Double> insideStromaKerObjectsDistance = []
        ArrayList<Double> outsideStromaKerObjectsDistance = []
        ArrayList<Double> insideStromaKerObjectsIntensity = []
        ArrayList<Double> outsideStromaKerObjectsIntensity = []

        selectObjectsByClassification("KER_488")
        def kerObjects = getSelectedObjects()
        kerObjects.each { obj ->
            def distance = measurement(obj, "Signed distance to annotation Stroma µm")
            def intensity = measurement(obj, "pNDRG1_647: Cell: Max")
            if (distance < 0) {
                insideStromaKerObjectsDistance.add(distance)
                insideStromaKerObjectsIntensity.add(intensity)
            } else if (distance > 0) {
                outsideStromaKerObjectsDistance.add(distance)
                outsideStromaKerObjectsIntensity.add(intensity)
            }
        }

        double[] insideStromaKerObjectsDistanceArray = insideStromaKerObjectsDistance as double[]
        double[] outsideStromaKerObjectsDistanceArray = outsideStromaKerObjectsDistance as double[]
        double[] insideStromaKerObjectsIntensityArray = insideStromaKerObjectsIntensity as double[]
        double[] outsideStromaKerObjectsIntensityArray = outsideStromaKerObjectsIntensity as double[]

        def pearsonInsideStromaKer = new PearsonsCorrelation().correlation(insideStromaKerObjectsDistanceArray, insideStromaKerObjectsIntensityArray)
        def pearsonOutsideStromaKer = new PearsonsCorrelation().correlation(outsideStromaKerObjectsDistanceArray, outsideStromaKerObjectsIntensityArray)

        println("KER_488 - Pearson correlation inside stroma: " + pearsonInsideStromaKer)
        println("KER_488 - Pearson correlation outside stroma: " + pearsonOutsideStromaKer)

        ArrayList<Double> insideStromaPNDRG1KerObjectsDistance = []
        ArrayList<Double> outsideStromaPNDRG1KerObjectsDistance = []
        ArrayList<Double> insideStromaPNDRG1KerObjectsIntensity = []
        ArrayList<Double> outsideStromaPNDRG1KerObjectsIntensity = []

        selectObjectsByClassification("KER_488: pNDRG1_647")
        def pNDRG1KerObjects = getSelectedObjects()
        pNDRG1KerObjects.each { obj ->
            def distance = measurement(obj, "Signed distance to annotation Stroma µm")
            def intensity = measurement(obj, "pNDRG1_647: Cell: Max")
            if (distance < 0) {
                insideStromaPNDRG1KerObjectsDistance.add(distance)
                insideStromaPNDRG1KerObjectsIntensity.add(intensity)
            } else if (distance > 0) {
                outsideStromaPNDRG1KerObjectsDistance.add(distance)
                outsideStromaPNDRG1KerObjectsIntensity.add(intensity)
            }
        }

        double[] insideStromaPNDRG1KerObjectsDistanceArray = insideStromaPNDRG1KerObjectsDistance as double[]
        double[] outsideStromaPNDRG1KerObjectsDistanceArray = outsideStromaPNDRG1KerObjectsDistance as double[]
        double[] insideStromaPNDRG1KerObjectsIntensityArray = insideStromaPNDRG1KerObjectsIntensity as double[]
        double[] outsideStromaPNDRG1KerObjectsIntensityArray = outsideStromaPNDRG1KerObjectsIntensity as double[]

        def pearsonInsideStromaKerPNDRG1 = new PearsonsCorrelation().correlation(insideStromaPNDRG1KerObjectsDistanceArray, insideStromaPNDRG1KerObjectsIntensityArray)
        def pearsonOutsideStromaKerPNDRG1 = new PearsonsCorrelation().correlation(outsideStromaPNDRG1KerObjectsDistanceArray, outsideStromaPNDRG1KerObjectsIntensityArray)

        println("KER_488: pNDRG1_647 - Pearson correlation inside stroma: " + pearsonInsideStromaKerPNDRG1)
        println("KER_488: pNDRG1_647 - Pearson correlation outside stroma: " + pearsonOutsideStromaKerPNDRG1)

        def folder_path = "/your/path/to/results"
        def imageName = getProjectEntry().getImageName()
        def path = buildFilePath(folder_path, "${imageName}.csv")

        if (!Files.exists(Paths.get(path))) {
            new PrintWriter(path).withWriter { writer ->
                writer.println("FN_568,sigma,pearsonInsideStromaKer,pearsonOutsideStromaKer,pearsonInsideStromaKerPNDRG1,pearsonOutsideStromaKerPNDRG1")
            }
        }

        new PrintWriter(new FileOutputStream(new File(path), true)).withWriter { writer ->
            writer.println("${fn_threshold},${sigma},${pearsonInsideStromaKer},${pearsonOutsideStromaKer},${pearsonInsideStromaKerPNDRG1},${pearsonOutsideStromaKerPNDRG1}")
        }

        println('All done!')
    }
}
