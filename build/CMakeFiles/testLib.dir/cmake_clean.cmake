file(REMOVE_RECURSE
  "libtestLib.pdb"
  "libtestLib.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/testLib.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
