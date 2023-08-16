// "Batch RGB Merge"

// Opens multiple sets of three separate color channels as
// an RGB stack or converts them to RGB images. File names
// ending in "d1", "d2" and "d0" are assumed to be red, green
// and blue channels respectively, but this can be changed in
// the dialog box.

// A sample image set, courtesy of Mikael Bjorklund, is available at:
//    http://rsb.info.nih.gov/ij/macros/images/DrosophilaCells.zip
// It consists of three images of Drosophila S2 cells,
// each with three channels (d0=blue, d1=red and d2=green  
// as indicated by the end of the filename. The staining is
// standard Hoechst, phalloidin, tubulin.

  Dialog.create("RGB Batch Convert");
  Dialog.addString("Green Suffix:", "green");
  Dialog.addString("Red Suffix:", "red");
  Dialog.addString("Blue Suffix:", "blue");
  Dialog.addString("Yellow Suffix:", "yellow");
  Dialog.addCheckbox("Open as Stack", true);
  Dialog.show();

  greenSuffix = Dialog.getString() + ".";
  redSuffix = Dialog.getString() + ".";
  blueSuffix = Dialog.getString() + ".";
  yellowSuffix = Dialog.getString() + ".";
  openAsStack = Dialog.getCheckbox();
  if (openAsStack)
      openImagesAsStack();
  else
      batchConvert();
  exit;

  function openImagesAsStack() {
      dir = getDirectory("Choose Source Directory ");
      list = getFileList(dir);
      setBatchMode(true);
      n = list.length;

      stack = 0;
      first = 0;
      for (i=0; i<n/4; i++) {
          showProgress(i+1, n/4);
          red="?"; green="?"; blue="?"; yellow="?";
          for (j=first; j<first+4; j++) {
              if (indexOf(list[j], redSuffix)!=-1)
                  red = list[j];
              if (indexOf(list[j], greenSuffix)!=-1)
                  green = list[j];
              if (indexOf(list[j], blueSuffix)!=-1)
                  blue = list[j];
              if (indexOf(list[j], yellowSuffix)!=-1)
                  yellow = list[j];
          }
          open(dir+red);
          open(dir+green);
          open(dir+blue);
          open(dir+yellow);
          run("RGB Merge...", "red=["+red+"] green=["+green+"] blue=["+blue+"] yellow=["+yellow+"]");
          width=getWidth; height=getHeight;
          run("Copy");
          close();
          if (stack==0) {
              newImage("RGB Stack", "RGB Black", width, height, n/4);
              stack = getImageID;
          }
          selectImage(stack);
          setSlice(i+1);
          run("Paste");
          index = indexOf(red, redSuffix);
          name = substring(red, 0, index);
          setMetadata(name);
          first += 4;
      }
      setSlice(1);
      run("Select None");
      setBatchMode(false);
  }

  function batchConvert() {
      dir1 = getDirectory("Choose Source Directory ");
      dir2 = getDirectory("Choose Destination Directory ");
      list = getFileList(dir1);
      setBatchMode(true);
      n = list.length;

      stack = 0;
      first = 0;
      for (i=0; i<n/4; i++) {
          showProgress(i+1, n/4);
          red="?"; green="?"; blue="?"; yellow="?";
          for (j=first; j<first+4; j++) {
              if (indexOf(list[j], redSuffix)!=-1)
                  red = list[j];
              if (indexOf(list[j], greenSuffix)!=-1)
                  green = list[j];
              if (indexOf(list[j], blueSuffix)!=-1)
                  blue = list[j];
              if (indexOf(list[j], yellowSuffix)!=-1)
                  yellow = list[j];
          }
          open(dir1 +red);
          open(dir1 +green);
          open(dir1 +blue);
          open(dir1 +yellow);
          run("RGB Merge...", "red=["+red+"] green=["+green+"] blue=["+blue+"] yellow=["+yellow+"]");
          index = indexOf(red, redSuffix);
          name = substring(red, 0, index);
          saveAs("tiff", dir2+name);
          first += 4;
      }
  }
