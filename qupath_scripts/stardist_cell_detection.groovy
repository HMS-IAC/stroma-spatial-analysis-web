/**
 * Authors: Pete Bankhead
 * Adapted by Antoine A. Ruzette, Simon F. Nørrelykke
 * Date: 2025-02-25
 *
 * This script detects cells using StarDist, 
 * and infers cell boundaries using an average cell radius.
 * It requires a StarDist frozen model (see .pb file) that takes a single channel as input.
 * 
 * Released under the MIT License (see LICENSE file)
 */

import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.gui.measure.ObservableMeasurementTableData;
import qupath.lib.scripting.QP
import qupath.ext.stardist.StarDist2D
import java.nio.file.Files
import java.nio.file.Paths
import groovy.time.*

setImageType('FLUORESCENCE');
setChannelNames(getCurrentImageData(),'DAPI','FITC KER','TRITC FN','CY5 pNDRG1')
// Select the user-created ROI
selectAnnotations();

// ---- STARDIST
// Detect all nuclei
// IMPORTANT! Replace this with the path to your StarDist model
// that takes a single channel as input (e.g. dsb2018_heavy_augment.pb)
// You can find some at https://github.com/qupath/models
// (Check credit & reuse info before downloading)

// modify the path to the model
def modelPath = "/Users/antoine/Harvard/IAC/image_analysis_paper_QP_projects/SUIT2_pNDRG1/project/classifiers/stardist/dsb2018_heavy_augment.pb"

//def avgNucleusRadius = 2.5
//def minNucleusArea = (Math.PI * (avgNucleusRadius)**2) / 1.5
def cellCytoplasmExpansion = 5 // SUIT2 cells
//def cellCytoplasmExpansion = 4 // AsPC cells

def getImagePixelSize() {
    def server = getCurrentServer()
    def cal = server.getMetadata().getPixelCalibration()
    def calibrationString = cal.toString()
    def xPixelSize = extractXPixelSize(calibrationString)
    println("Pixel size [microns]: ${xPixelSize}")
    return xPixelSize as Double // Explicitly convert to double
}

def extractXPixelSize(calibrationString) {
    def pattern = /x=(\d+(\.\d+)?)\sµm/
    def matcher = (calibrationString =~ pattern)

    if (matcher.find()) {
        return matcher.group(1).toDouble()
    } else {
        return null
    }
}

def xPixelSize = getImagePixelSize()
def cellExpansionPixels = Math.round(cellCytoplasmExpansion / xPixelSize)

def dnn = DnnTools.builder(modelPath).build();
def stardist = StarDist2D.builder(dnn)
    .channels('DAPI')  
    .threshold(0.6)              // Prediction threshold
    .normalizePercentiles(1,99) // Percentile normalization. REQUIRED FOR IMC DATA
    .pixelSize(xPixelSize)              // Resolution for detection
    .includeProbability(true)
    .measureIntensity()
    .measureShape()
    .cellExpansion(cellExpansionPixels) // in pixels; 
    .constrainToParent(false)
    .build()

// Define which objects will be used as the 'parents' for detection
// Use QP.getAnnotationObjects() if you want to use all annotations, rather than selected objects
def pathObjects = QP.getSelectedObjects()

// Run detection for the selected objects
def imageData = QP.getCurrentImageData()
if (pathObjects.isEmpty()) {
    QP.getLogger().error("No parent objects are selected!")
    return
}

// Get cell detections
def cells = getCellObjects()

// Check if cells have already been detected
if (cells.isEmpty()) {
    // QP.getLogger().error("StarDist has not ran yet!")
    println('Starting cell detection...')
    def timeStart_CellDetection = new Date()
    stardist.detectObjects(imageData, pathObjects)
    dnn.getPredictionFunction().net.close()
    stardist.close() // This can help clean up & regain memory
    TimeDuration CellDetection_duration = TimeCategory.minus(new Date(), timeStart_CellDetection)
    println("Stardist running time: ${CellDetection_duration}")

    // Remove nuclei with extremely small areas
    def detections = getDetectionObjects()
    def cellAreas = detections.collect { measurement(it, "Nucleus: Area µm^2") }
    // Sort the total intensities
    Collections.sort(cellAreas)
    // Calculate the 5th and 95th percentiles
    def areaIndex5th = Math.round(cellAreas.size() * 0.05).toInteger()
    def areaPercentile5th = cellAreas[areaIndex5th]
    def toDeleteBelow5th = detections.findAll {measurement(it, 'Nucleus: Area µm^2') < areaPercentile5th}
    removeObjects(toDeleteBelow5th, true)

    println('Done!')
} else {
    println('StarDist already ran. Using existing cell detections!')
}