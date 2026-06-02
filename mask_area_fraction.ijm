// =====================================================
// Generic FIJI macro: measure mask area in batch
// =====================================================

// ---------------- CONFIG ----------------
inputPath = "/path/to/your/data/";
outputFile = File.getParent(inputPath) + File.separator + "area_analysis.csv";

fileSuffix = "_mask.tif";

useCircularROI = "NO";   // "YES" or "NO" - do you want to only analyze a centered circle? (for full well captures)
roiDiameter = 6000; // diameter of the centered circle if "YES" above

// Measurement settings
decimalPlaces = 3;
// ----------------------------------------

run("Set Measurements...", "area area_fraction limit redirect=None decimal=" + decimalPlaces);

// Overwrite output CSV and write header
if (File.exists(outputFile)) File.delete(outputFile);
File.append("Filename,Area,%Area\n", outputFile);

fileList = getFileList(inputPath);

for (i = 0; i < fileList.length; i++) {
    if (endsWith(fileList[i], fileSuffix)) {
        measureMask(inputPath + fileList[i], fileList[i]);
    }
}

print("Saved CSV to: " + outputFile);
print("All done!");


// =====================================================
// Measure one file
// =====================================================
function measureMask(fullPath, filename) {
    open(fullPath);
    print("Processing: " + filename);

    width = getWidth();
    height = getHeight();

    if (useCircularROI == "YES") {
        radius = roiDiameter / 2;
        x = width / 2 - radius;
        y = height / 2 - radius;
        makeOval(x, y, roiDiameter, roiDiameter);
    } else {
        // Measure whole image
        makeRectangle(0, 0, width, height);
    }

    run("Measure");

    area = getResult("Area", nResults - 1);
    areaFrac = getResult("%Area", nResults - 1);

    line = filename + "," + d2s(area, decimalPlaces) + "," + d2s(areaFrac, decimalPlaces);
    File.append(line, outputFile);

    run("Clear Results");
    close();
}