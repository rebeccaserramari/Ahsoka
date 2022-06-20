add_library(AhsokaLib SHARED 
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

target_link_libraries(Ahsoka AhsokaLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(Ahsoka AhsokaLib ${JELLYFISH_LIBRARIES})


#target_link_libraries(Ahsoka AhsokaLib ${WHATSHAP_LDFLAGS_OTHER})
#target_link_libraries(Ahsoka AhsokaLib ${WHATSHAP_LIBRARIES})
