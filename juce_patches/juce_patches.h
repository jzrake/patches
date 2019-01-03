/******************************************************************************
BEGIN_JUCE_MODULE_DECLARATION

  ID:               juce_patches
  vendor:           Jonathan Zrake
  version:          0.0.1
  name:             patches
  description:      A C++14 micro-library for adaptive mesh refinement (AMR) in computational gasdynamics applications.
  website:          https://github.com/jzrake/patches
  license:          GPL
  dependencies:     juce_ndarray

END_JUCE_MODULE_DECLARATION
******************************************************************************/


#pragma once
#define JUCE_PATCHES_INCLUDED
#include <juce_ndarray/juce_ndarray.h>
#include <juce_patches/src/patches.hpp>
#include <juce_patches/src/serializer.hpp>
