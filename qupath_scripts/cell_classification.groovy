/**
 * Authors: Antoine A. Ruzette, Simon F. Nørrelykke
 * Date: 2024-03-03
 *
 * This scripts classifies cells based on intensity thresholds of two markers.
 * The script is designed to be run on detections.
 * 
 * Released under the MIT License (see LICENSE file)
 */

import qupath.lib.objects.PathObject
import static qupath.lib.gui.scripting.QPEx.*
import java.nio.file.Files
import java.nio.file.Paths

// Iterate over threshold combinations
//def ker_thresholds = (0..2000000).step(100000).collect()
//def marker_thresholds = (0..500000).step(25000).collect()

 def ker_thresholds = [1200]
 def marker_thresholds = [650]

// Get the objects (here, we use detections - change if required)
def pathObjects = getDetectionObjects()

// Loop over each combination of thresholds
ker_thresholds.each { ker_threshold ->
    marker_thresholds.each { marker_threshold ->
        println('Starting cell sub-classification...')
        println('Resetting detections class...')
        resetDetectionClassifications()
        println('Parameters: ')
        print("Keratin threshold: ${ker_threshold}; Marker threshold: ${marker_threshold}")

        // Define cell classifier JSON string with updated thresholds
        def cellClassifierString = """
        {
          "object_classifier_type": "CompositeClassifier",
          "classifiers": [
            {
              "object_classifier_type": "SimpleClassifier",
              "function": {
                "classifier_fun": "ClassifyByMeasurementFunction",
                "measurement": "FITC KER: Cytoplasm: Median",
                "pathClassEquals": "FITC KER",
                "pathClassAbove": "FITC KER",
                "threshold": ${ker_threshold}
              },
              "pathClasses": [
                "FITC KER"
              ],
              "filter": "CELLS",
              "timestamp": 1694426186386
            },
            {
              "object_classifier_type": "SimpleClassifier",
              "function": {
                "classifier_fun": "ClassifyByMeasurementFunction",
                "measurement": "CY5 pNDRG1: Cell: Max",
                "pathClassEquals": "CY5 pNDRG1",
                "pathClassAbove": "CY5 pNDRG1",
                "threshold": ${marker_threshold}
              },
              "pathClasses": [
                "CY5 pNDRG1"
              ],
              "filter": "CELLS",
              "timestamp": 1694426219269
            }
          ]
        }
        """

        // Define the file path to save the JSON, including the value of thresholds
        def cellClassifierFilePath = "./classifiers/cell_classifiers/cell_classifier_ker-thresh${ker_threshold}_marker-thresh${marker_threshold}.json"

        // Save the JSON string to a JSON file
        def cellClassifierFile = new File(cellClassifierFilePath)
        cellClassifierFile.text = cellClassifierString

        // Execute your classifier with the current thresholds
        runObjectClassifier(cellClassifierFilePath)

        // Update hierarchy
        fireHierarchyUpdate()
    }
}