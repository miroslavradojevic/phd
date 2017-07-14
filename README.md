# Neuron tracing using probability hypothesis filtering #

ImageJ plugin, [https://doi.org/10.1093/bioinformatics/btw751](Link URL)

### Usage ###

* copy [phd_.jar](https://bitbucket.org/miroslavradojevic/phd/src) to ImageJ plugins directory
* copy [vess_.jar](https://bitbucket.org/miroslavradojevic/vess/src) to ImageJ plugins directory
* Open ImageJ > BrainCadet > PHD
* select file (works with 8bit images and image stacks)
* choose values in the parameter menu

### ImageJ macro for batch processing ###

* parameter grid 
* possibility to use the parameter grid to process all existing *.tif 8bit images in given directory 
* ij macro: ij-macro/run_phd.ijm
* test case (~/test/ directory with run_phd.ijm and img.tif), example call: 
```
#!java

java -Xmx4g -jar ~/ImageJ/ij.jar -ijpath ~/ImageJ/plugins/ -batch ~/test/run_phd.ijm ~/test/
```