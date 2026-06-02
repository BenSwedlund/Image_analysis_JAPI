// =====================================================
// Generic FIJI macro: create binary masks from multi-channel fluorescent images
// =====================================================

// ---------------- CONFIG ----------------
path = "/path/to/your/data/";
batch_mode = "NO";           // "YES" or "NO"
inputImage = "test.tif";
tomatch = "";                // optional filename filter

channelsToProcess = newArray("1", "2"); //exclude brightfield in 3

// Channel labels used in output filenames
channelNames = newArray("mCherry", "GFP");

// Per-channel threshold ranges
min_t_arr  = newArray(250, 750);
max_t_arr  = newArray(65535, 65535);

min_t2_arr  = newArray(100, 100);
max_t2_arr  = newArray(255, 255);

// Optional display stretch for channel previews
set_min_arr = newArray(100, 100);
set_max_arr = newArray(500, 2000);

// Crop settings - do you want to crop? If yes, is it a centered or off-centered square?
// crop_x and crop_y are only used if the crop is not centered
use_center_crop = "YES";
crop_diameter = 2500;

crop_x = 1060;
crop_y = 990;


// Blurring and eroding parameters - important as it sets the graularity of the mask
dilate_erode_iterations = 3;
first_gaussian_blur = 8;
second_gaussian_blur = 8;

// Threshold behavior - automatic or set?
use_manual_first_threshold = "YES";
use_manual_second_threshold = "YES";
invert_LUT = "YES";

// ----------------------------------------

if (batch_mode == "YES") {
    fileList = getFileList(path);
    for (f = 0; f < fileList.length; f++) {
        if ((endsWith(fileList[f], ".tif") || endsWith(fileList[f], ".tiff")) &&
            (tomatch == "" || indexOf(fileList[f], tomatch) != -1)) {
            processFile(fileList[f]);
        }
    }
} else {
    processFile(inputImage);
}

print("All done!");


// =====================================================
// Main processing
// =====================================================
function processFile(fileName) {
    fullPath = path + fileName;
    baseName = getBaseName(fileName);

    open(fullPath);
    origTitle = getTitle();

    if (use_center_crop == "YES") {
        width = getWidth();
        height = getHeight();
        crop_x = (width - crop_diameter) / 2;
        crop_y = (height - crop_diameter) / 2;
    }
    makeRectangle(crop_x, crop_y, crop_diameter, crop_diameter);

    outPath = path + "mask/";
    outPath_masks = outPath + "TIF-masks/";
    makeDir(outPath);
    makeDir(outPath_masks);

    for (c = 0; c < lengthOf(channelsToProcess); c++) {
        channel = channelsToProcess[c];
        channelName = channelNames[c];

        processChannel(origTitle, baseName, channel, channelName, c, outPath, outPath_masks);
    }

    run("Close All");
}


// =====================================================
// Per-channel processing
// =====================================================
function processChannel(origTitle, baseName, channel, channelName, c, outPath, outPath_masks) {
    imageName = baseName + "_" + channelName + ".tif";
    maskName = baseName + "_" + channelName + "_mask.tif";

    print("Processing channel " + channel + " -> " + channelName);

    selectWindow(origTitle);
    run("Duplicate...", "title=" + imageName + " duplicate channels=" + channel);

    selectWindow(imageName);
    run("Duplicate...", "title=" + maskName);

    // --- Make binary mask ---
    selectWindow(maskName);
    setOption("BlackBackground", false);

    if (use_manual_first_threshold == "YES") {
        setThreshold(min_t_arr[c], max_t_arr[c]);
    } else {
        setAutoThreshold("Default");
    }
    run("Convert to Mask");

    run("Gaussian Blur...", "sigma=" + first_gaussian_blur);

    if (use_manual_second_threshold == "YES") {
        setThreshold(min_t2_arr[c], max_t2_arr[c]);
    } else {
        setAutoThreshold("Default");
    }
    run("Convert to Mask");

    run("Gaussian Blur...", "sigma=" + second_gaussian_blur);
    setAutoThreshold("Default");
    run("Convert to Mask");

    for (i = 0; i < dilate_erode_iterations; i++) run("Dilate");
    for (i = 0; i < dilate_erode_iterations; i++) run("Erode");

    saveAs("Tiff", outPath_masks + maskName);

    if (invert_LUT == "YES") run("Invert LUT");
    saveAs("Jpeg", outPath + baseName + "_" + channelName + "_mask.jpeg");

    // --- Save channel preview ---
    selectWindow(imageName);
    setMinAndMax(set_min_arr[c], set_max_arr[c]);
    saveAs("Jpeg", outPath + baseName + "_" + channelName + ".jpeg");

    // --- Overlay mask on channel preview ---
    run("Add Image...", "image=" + maskName + " x=0 y=0 opacity=30");
    saveAs("Jpeg", outPath + baseName + "_" + channelName + "_withmask.jpeg");
}


// =====================================================
// Helpers
// =====================================================
function makeDir(p) {
    if (!File.exists(p)) File.makeDirectory(p);
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