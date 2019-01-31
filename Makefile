all: juce_patches/src/patches.cpp juce_patches/src/patches.hpp

juce_patches/src/patches.cpp: patches.cpp
	cp $^ $@

juce_patches/src/patches.hpp: patches.hpp
	sed '/#include "ndarray.hpp"/d' $^ > $@
