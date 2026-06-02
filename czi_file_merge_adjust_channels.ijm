// ============================================================
// Generic CZI batch processing template
// Edit only the CONFIG section for new datasets
// ============================================================

// -------------------- CONFIG --------------------
basePath = "/path/to/your/data/";
folderNames = newArray("Day4");          // Name of subfolders to analyze, e.g. newArray("Day1", "Day2")
inputExtension = ".czi";

// Channel mapping: order in merged image
channelOrder = newArray("C=2", "C=0", "C=1");
channelColors = newArray("c1", "c2", "c4"); // display colors in Merge Channels (c1 = red, c2 = green, c3 = blue, c4 = gray)

// Output labels of the different channels (e.g. GFP, mCherry, brightfield)
outputLabels = newArray("mCherry", "GFP", "brightfield");

// Background subtraction, size of Rolling ball
backgroundRolling = 50;

// Display scaling
lowerPercentile1 = 40;
upperPercentile1 = 95;
lowerPercentile2 = 30;
upperPercentile2 = 95;

// Thresholds for deciding whether to use fixed high max
threshPos1 = 2000;
pctSwitch1 = 0.1;

threshPos2 = 200;
pctSwitch2 = 0.1;

// Brightfield settings / third channel display range
minCh3 = 0;
maxCh3 = 16500;

// Output format names
saveMergedTif = true;
saveAllJpeg = true;
saveSingleChannelJpegs = true;

// ------------------------------------------------

for (f = 0; f < folderNames.length; f++) {
    processFolder(basePath, folderNames[f]);
}

print("All files processed.");


// ============================================================
// Folder processing
// ============================================================
function processFolder(basePath, folderName) {
    inputDir = basePath + folderName + "/";
    outputBase = inputDir + "Processed/";
    outputTif = outputBase + "TIF/";
    outputJpeg = outputBase + "JPEG/";
    outputJpegAll = outputJpeg + "all/";

    makeDir(outputBase);
    makeDir(outputTif);
    makeDir(outputJpeg);
    makeDir(outputJpegAll);

    list = getFileList(inputDir);

    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], inputExtension)) {
            processCziImage(inputDir + list[i], outputTif, outputJpeg, outputJpegAll);
        }
    }
}


// ============================================================
// Image processing
// ============================================================
function processCziImage(inputImage, outputTif, outputJpeg, outputJpegAll) {
    baseName = getBaseName(inputImage);

    options = "autoscale color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT series=0";
    run("Bio-Formats Importer", "open=[" + inputImage + "] " + options);

    // Background subtraction on selected channels
    channelsToCorrect = Array.slice(channelOrder, 0, channelOrder.length - 1);
    for (i = 0; i < channelsToCorrect.length; i++) {
        title = baseName + ".czi - " + baseName + ".czi #1 - " + channelsToCorrect[i];
        selectImage(title);
        run("Subtract Background...", "rolling=" + backgroundRolling);
    }

    // Merge channels
    ch1 = baseName + ".czi - " + baseName + ".czi #1 - " + channelOrder[0];
    ch2 = baseName + ".czi - " + baseName + ".czi #1 - " + channelOrder[1];
    ch3 = baseName + ".czi - " + baseName + ".czi #1 - " + channelOrder[2];

    run("Merge Channels...", channelColors[0] + "=[" + ch1 + "] " +
                            channelColors[1] + "=[" + ch2 + "] " +
                            channelColors[2] + "=[" + ch3 + "] create");

    // Scale channel 1
    setSlice(1);
    frac1 = getSignalFraction(threshPos1);
    print(baseName + " signal1 fraction = " + frac1);

    if (frac1 < pctSwitch1) {
        upper = getPercentileValue(99);
    } else {
        upper = getPercentileValue(upperPercentile1);
    }
    lower = getPercentileValue(lowerPercentile1);
    setMinAndMax(lower, upper);

    // Scale channel 2
    setSlice(2);
    frac2 = getSignalFraction(threshPos2);
    print(baseName + " signal2 fraction = " + frac2);

    if (frac2 < pctSwitch2) {
        upper = getPercentileValue(99);
    } else {
        upper = getPercentileValue(upperPercentile2);
    }
    lower = getPercentileValue(lowerPercentile2);
    setMinAndMax(lower, upper);

    // Scale channel 3
    setSlice(3);
    setMinAndMax(minCh3, maxCh3);

    // Add scale bar
    run("Scale Bar...", "width=100 height=100 horizontal bold overlay");

    // Save outputs
    if (saveMergedTif)
        saveAs("Tiff", outputTif + baseName + "_all.tif");

    if (saveAllJpeg)
        saveAs("Jpeg", outputJpegAll + baseName + "_all.jpeg");

    if (saveSingleChannelJpegs) {
        Stack.setActiveChannels("100");
        saveAs("Jpeg", outputJpeg + baseName + "_signal1.jpg");

        Stack.setActiveChannels("010");
        saveAs("Jpeg", outputJpeg + baseName + "_signal2.jpg");
    }

    close();
}


// ============================================================
// Helper functions
// ============================================================
function makeDir(path) {
    if (!File.exists(path)) File.makeDirectory(path);
}

function getBaseName(filePath) {
    parts = split(filePath, File.separator);
    title = parts[parts.length - 1];
    dotIndex = lastIndexOf(title, ".");
    if (dotIndex > 0)
        return substring(title, 0, dotIndex);
    else
        return title;
}

function getPercentileValue(p) {
    getRawStatistics(nPixels, mean, min, max, std, histogram);

    nBins = histogram.length;
    total = 0;
    for (i = 0; i < nBins; i++) total += histogram[i];

    target = total * (p / 100.0);
    cumulative = 0;

    for (i = 0; i < nBins; i++) {
        cumulative += histogram[i];
        if (cumulative >= target) {
            return min + (i / (nBins - 1)) * (max - min);
        }
    }
    return max;
}

function getSignalFraction(threshPos) {
    getRawStatistics(nPixels, mean, min, max, std, histogram);
    nBins = histogram.length;

    idx = round((threshPos - min) * (nBins - 1) / (max - min));
    if (idx < 0) idx = 0;
    if (idx >= nBins) idx = nBins - 1;

    above = 0;
    for (i = idx; i < nBins; i++)
        above += histogram[i];

    return above / nPixels;
}