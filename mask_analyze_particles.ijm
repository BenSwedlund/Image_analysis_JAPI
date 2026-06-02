// =====================================================
// Generic FIJI macro: particle analysis on mask images
// =====================================================

// ---------------- CONFIG ----------------
dir = "/path/to/your/data/";
fileSuffix = "_all_mask.tif";

// Optional ROI behavior
useCircularROI = "NO";   // "YES" or "NO"
roiDiameter = 6000;      // used only if useCircularROI = "YES"

// Output folder for particle tables
outDir = dir + "particles/";
File.makeDirectory(outDir);

// Measurement settings
measurementOptions = "area centroid fit limit redirect=None decimal=3";

// Particle analysis settings
particleOptions = "show=Overlay display include add";
// ----------------------------------------

// Get files
list = getFileList(dir);

for (i = 0; i < list.length; i++) {
    filename = list[i];

    if (endsWith(filename, fileSuffix)) {
        processMask(dir + filename, filename, outDir, measurementOptions, particleOptions);
    }
}

print("All done!");


// =====================================================
// Process one mask file
// =====================================================
function processMask(fullPath, filename, outDir, measurementOptions, particleOptions) {
    open(fullPath);
    print("Processing: " + filename);

    // Clear results before each image
    run("Set Measurements...", measurementOptions);
    run("Clear Results");

    width = getWidth();
    height = getHeight();

    if (useCircularROI == "YES") {
        radius = roiDiameter / 2;
        x = width / 2 - radius;
        y = height / 2 - radius;
        makeOval(x, y, roiDiameter, roiDiameter);
    } else {
        makeRectangle(0, 0, width, height);
    }

    // Run particle analysis
    run("Analyze Particles...", particleOptions);

    // Add filename column to every row in the current Results table
    n = nResults;
    for (j = 0; j < n; j++) {
        setResult("Filename", j, filename);
    }
    updateResults();

    // Save results
    outCsv = outDir + replace(filename, ".tif", "_particles.csv");
    saveAs("Results", outCsv);

    close();
    run("Clear Results");
}