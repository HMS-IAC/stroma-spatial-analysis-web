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
def sigmas = 15
def fn_thresholds = 5629.18

//def sigmas = [25]
//def fn_thresholds = [3500]

// ---- STROMA ANNOTATION
println('Starting stroma annotation...')

// Loop over each combination of thresholds
sigmas.each { sigma ->
    fn_thresholds.each { fn_threshold ->


        // ---- REMOVE STROMA ANNOTATION FOR NEXT ITERATION
        
        // Define the name of the annotation to be removed
        selectObjects { p -> p.getPathClass() == getPathClass("Stroma") && p.isAnnotation() }
        
        def selectedObjects = QP.getSelectedObjects()
        
        if (!selectedObjects.isEmpty()) {
            // Remove the selected annotations
            removeObjects(selectedObjects, false)
            println("Selected Stroma annotations deleted.")
        } else {
            println("No Stroma annotations selected for deletion.")
        }
    
        resetSelection();
        // Select the user-created ROI
        selectAnnotations();
        
        println('Parameters:')
        println("Sigma: ${sigma}; FN_568_threshold: ${fn_threshold}")
        
        // Create a JSON string with sigmaX included
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
        // Define the file path to save the JSON, including the value of sigmaX
        def pixelClassifierFilePath = "/Users/antoine/Harvard/IAC/image_analysis_paper_QP_projects/SUIT2_pNDRG1/project/thresholds/classifiers/stroma_annotation/stroma_annotator_sigma=${sigma}_fn-thresh=${fn_threshold}.json"
        
        // Save the JSON string to a JSON file
        def pixelClassifierFile = new File(pixelClassifierFilePath)
        pixelClassifierFile.text = pixelClassifierString
        
        createAnnotationsFromPixelClassifier(pixelClassifierFilePath, 0.0, 0.0)
        println('Done!')
        
        // ---- SIGNED DISTANCE CALCULATION
        
        println('Starting calculation of signed distance...')
        detectionToAnnotationDistancesSigned(false)
        println('Done!')


        // ---- CALCULATE CORRELATION BEFORE AND AFTER DISTANCE = 0
        // ---- KER_488

        println("Starting calculation of Pearson's correlation coefficients...")

        // Initialize lists to store objects
        ArrayList<Double> insideStromaKerObjectsDistance = []
        ArrayList<Double> outsideStromaKerObjectsDistance = []
        ArrayList<Double> insideStromaKerObjectsIntensity = []
        ArrayList<Double> outsideStromaKerObjectsIntensity = []

        // Select KER_488 objects and iterate over them
        selectObjectsByClassification("KER_488")
        kerObjects = getSelectedObjects()
        kerObjects.each { obj ->
            def distance = measurement(obj, "Signed distance to annotation Stroma µm")
            def intensity = measurement(obj, "pNDRG1_647: Cell: MedianTotalIntensity")
            if (distance < 0) {
                insideStromaKerObjectsDistance.add(distance)
                insideStromaKerObjectsIntensity.add(intensity)
            } else if (distance > 0) {
                outsideStromaKerObjectsDistance.add(distance)
                outsideStromaKerObjectsIntensity.add(intensity)
            }
        }

//        // Convert ArrayList to arrays
//        double[] insideStromaKerObjectsDistanceArray = insideStromaKerObjectsDistance.stream().mapToDouble{ it }.toArray()
//        double[] outsideStromaKerObjectsDistanceArray = outsideStromaKerObjectsDistance.stream().mapToDouble{ it }.toArray()
//        double[] insideStromaKerObjectsIntensityArray = insideStromaKerObjectsIntensity.stream().mapToDouble{ it }.toArray()
//        double[] outsideStromaKerObjectsIntensityArray = outsideStromaKerObjectsIntensity.stream().mapToDouble{ it }.toArray()
//
//        // Calculate Pearson correlation for inside stroma
//        def pearsonInsideStromaKer = new PearsonsCorrelation().correlation(insideStromaKerObjectsDistanceArray, insideStromaKerObjectsIntensityArray)
//        def pearsonOutsideStromaKer = new PearsonsCorrelation().correlation(outsideStromaKerObjectsDistanceArray, outsideStromaKerObjectsIntensityArray)
//        // Output correlations
//        println("KER_488 - Pearson correlation inside stroma: " + pearsonInsideStromaKer)
//        println("KER_488 - Pearson correlation outside stroma: " + pearsonOutsideStromaKer)
//
//
//        // ---- CALCULATE CORRELATION BEFORE AND AFTER DISTANCE = 0
//        // ---- KER_488: pNDRG1_647
//        
//        ArrayList<Double> insideStromaPNDRG1KerObjectsDistance = []
//        ArrayList<Double> outsideStromaPNDRG1KerObjectsDistance = []
//        ArrayList<Double> insideStromaPNDRG1KerObjectsIntensity = []
//        ArrayList<Double> outsideStromaPNDRG1KerObjectsIntensity = []
//
//        // Select KER_488: pNDRG1_647 objects and iterate over them
//        selectObjectsByClassification("KER_488: pNDRG1_647")
//        pNDRG1KerObjects = getSelectedObjects()
//        pNDRG1KerObjects.each { obj ->
//            def distance = measurement(obj, "Signed distance to annotation Stroma µm")
//            def intensity = measurement(obj, "pNDRG1_647: Cell: MedianTotalIntensity")
//            if (distance < 0) {
//                insideStromaPNDRG1KerObjectsDistance.add(distance)
//                insideStromaPNDRG1KerObjectsIntensity.add(intensity)
//            } else if (distance > 0) {
//                outsideStromaPNDRG1KerObjectsDistance.add(distance)
//                outsideStromaPNDRG1KerObjectsIntensity.add(intensity)
//            }
//        }
//        
//        // Convert lists to arrays
//        double[] insideStromaPNDRG1KerObjectsDistanceArray = insideStromaPNDRG1KerObjectsDistance as double[]
//        double[] outsideStromaPNDRG1KerObjectsDistanceArray = outsideStromaPNDRG1KerObjectsDistance as double[]
//        double[] insideStromaPNDRG1KerObjectsIntensityArray = insideStromaPNDRG1KerObjectsIntensity as double[]
//        double[] outsideStromaPNDRG1KerObjectsIntensityArray = outsideStromaPNDRG1KerObjectsIntensity as double[]
//
//        // Calculate Pearson correlation for inside stroma
//        def pearsonInsideStromaKerPNDRG1 = new PearsonsCorrelation().correlation(insideStromaPNDRG1KerObjectsDistanceArray, insideStromaPNDRG1KerObjectsIntensityArray)
//        def pearsonOutsideStromaKerPNDRG1 = new PearsonsCorrelation().correlation(outsideStromaPNDRG1KerObjectsDistanceArray, outsideStromaPNDRG1KerObjectsIntensityArray)
//
//        // Calculate Spearman correlation
//        def spearmanInsideStromaKer = new SpearmansCorrelation().correlation(insideStromaKerObjectsDistanceArray, insideStromaKerObjectsIntensityArray)
//        def spearmanOutsideStromaKer = new SpearmansCorrelation().correlation(outsideStromaKerObjectsDistanceArray, outsideStromaKerObjectsIntensityArray)
//        def spearmanInsideStromaKerPNDRG1 = new SpearmansCorrelation().correlation(insideStromaPNDRG1KerObjectsDistanceArray, insideStromaPNDRG1KerObjectsIntensityArray)
//        def spearmanOutsideStromaKerPNDRG1 = new SpearmansCorrelation().correlation(outsideStromaPNDRG1KerObjectsDistanceArray, outsideStromaPNDRG1KerObjectsIntensityArray)
//
//        // Output correlations
//        println("KER_488: pNDRG1_647 - Pearson correlation inside stroma: " + pearsonInsideStromaKerPNDRG1)
//        println("KER_488: pNDRG1_647 - Pearson correlation outside stroma: " + pearsonOutsideStromaKerPNDRG1)
//
//        println('Done!')
//
//        // Export Class Ratio and corresponding parameters
//        println('Starting exportation of parameters and correlations to CSV...')
//        def folder_path = "/Users/antoine/Harvard/IAC/NinaKozlova/20240111_AsPC_pNDRG1/Results/20240229_correlations"
//        def imageName = getProjectEntry().getImageName()
//
//        def path = buildFilePath(folder_path, "${imageName}" + ".csv")
//        def measurements = ["FN_568", "sigma", "pearsonInsideStromaKer", "pearsonOutsideStromaKer", "pearsonInsideStromaKerPNDRG1", "pearsonOutsideStromaKerPNDRG1", "spearmanInsideStromaKer", "spearmanOutsideStromaKer", "spearmanInsideStromaKerPNDRG1", "spearmanOutsideStromaKerPNDRG1"]
//
//        // Check if the file exists, if not, create a new file with a header
//        if (!Files.exists(Paths.get(path))) {
//            try (def writer = new PrintWriter(path)) {
//                def sb = new StringBuilder()
//                sb.append("FN_568, sigma, pearsonInsideStromaKer, pearsonOutsideStromaKer, pearsonInsideStromaKerPNDRG1, pearsonOutsideStromaKerPNDRG1, spearmanInsideStromaKer, spearmanOutsideStromaKer, spearmanInsideStromaKerPNDRG1, spearmanOutsideStromaKerPNDRG1")
//                writer.println(sb.toString())
//            }
//        }
//
//        // Append data to the CSV file
//        try (def writer = new PrintWriter(new FileOutputStream(new File(path), true))) {
//            def sb = new StringBuilder()
//            sb.append(fn_threshold)
//            sb.append(',')
//            sb.append(sigma)
//            sb.append(',')
//            sb.append(pearsonInsideStromaKer)
//            sb.append(',')
//            sb.append(pearsonOutsideStromaKer)
//            sb.append(',')
//            sb.append(pearsonInsideStromaKerPNDRG1)
//            sb.append(',')
//            sb.append(pearsonOutsideStromaKerPNDRG1)
//            sb.append(',')
//            sb.append(spearmanInsideStromaKer)
//            sb.append(',')
//            sb.append(spearmanOutsideStromaKer)
//            sb.append(',')
//            sb.append(spearmanInsideStromaKerPNDRG1)
//            sb.append(',')
//            sb.append(spearmanOutsideStromaKerPNDRG1)
//            writer.println(sb.toString())
//        } catch (Exception e) {
//            println "An error occurred while writing to the CSV file: ${e.message}"
//        }
//        println('All done!')

    }
}
        