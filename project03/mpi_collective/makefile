build:
	rm -f cs
	mpic++ -O2 -Wall main.cpp cs.cpp bitmap.cpp debug.cpp -Wno-unused-but-set-variable -Wno-unused-function -Wno-write-strings -Wno-unused-result -o cs

run:
	mpiexec -n 8 cs sunset.bmp temp.bmp 13
	diff temp.bmp stretched-13.bmp