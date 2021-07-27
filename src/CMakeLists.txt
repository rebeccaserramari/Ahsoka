add_library(PolyassemblyLib SHARED 
				polyassembly.cpp
				argumentparser.cpp
				timer.cpp
				graph.cpp
				alignmentreader.cpp
				jellyfishkmercounter.cpp
				uniquekmers.cpp
				
				/home/rebecca/work/whatshap-code/whatshap/src/read.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/readset.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/entry.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/indexset.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/polyphase/readscoring.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/polyphase/trianglesparsematrix.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/polyphase/clustereditingsolution.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/polyphase/clustereditingsolver.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/polyphase/inducedcostheuristic.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/polyphase/staticsparsegraph.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/polyphase/edgeheap.cpp
				/home/rebecca/work/whatshap-code/whatshap/src/polyphase/haplothreader.cpp
				)
				
				
add_executable(Ahsoka polyassembly.cpp)

target_link_libraries(Ahsoka PolyassemblyLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(Ahsoka PolyassemblyLib ${JELLYFISH_LIBRARIES})

#target_link_libraries(Ahsoka PolyassemblyLib ${WHATSHAP_LDFLAGS_OTHER})
#target_link_libraries(Ahsoka PolyassemblyLib ${WHATSHAP_LIBRARIES})