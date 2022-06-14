add_library(PolyassemblyLib SHARED 
				polyassembly.cpp
				argumentparser.cpp
				graph.cpp
				alignmentreader.cpp
				
				
				~/whatshap-code/src/read.cpp
				~/whatshap-code/src/readset.cpp
				~/whatshap-code/src/entry.cpp
				~/whatshap-code/src/indexset.cpp
				~/whatshap-code/src/polyphase/readscoring.cpp
				~/whatshap-code/src/polyphase/trianglesparsematrix.cpp
				~/whatshap-code/src/polyphase/clustereditingsolution.cpp
				~/whatshap-code/src/polyphase/clustereditingsolver.cpp
				~/whatshap-code/src/polyphase/inducedcostheuristic.cpp
				~/whatshap-code/src/polyphase/staticsparsegraph.cpp
				~/whatshap-code/src/polyphase/edgeheap.cpp
				~/whatshap-code/src/polyphase/haplothreader.cpp
				)
				
				
add_executable(Ahsoka polyassembly.cpp)

target_include_directories(PolyassemblyLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Ahsoka PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(Ahsoka PolyassemblyLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(Ahsoka PolyassemblyLib ${JELLYFISH_LIBRARIES})

#target_link_libraries(Ahsoka PolyassemblyLib ${WHATSHAP_LDFLAGS_OTHER})
#target_link_libraries(Ahsoka PolyassemblyLib ${WHATSHAP_LIBRARIES})
