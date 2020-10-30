file(REMOVE_RECURSE
  "liblibC.pdb"
  "liblibC.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/libC.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
