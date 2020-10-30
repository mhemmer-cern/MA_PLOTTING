file(REMOVE_RECURSE
  "liblibB.pdb"
  "liblibB.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/libB.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
