all: display-image

display-image: display-image.cc
	gcc -v -std=c++11 -x c++ display-image.cc -I/usr/local/include/ -L/usr/lib/ -lstdc++ -L/usr/local/lib -lopencv_highgui -lopencv_core -lopencv_imgcodecs -o display-image

clean:
	rm -v display-image


