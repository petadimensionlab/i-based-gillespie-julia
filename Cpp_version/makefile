MAIN_CPP = blowfly
ifeq ($(OS),Windows_NT)
    CXX = g++
    TARGET = a.exe
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        #CXX = /usr/local/bin/g++
        #CXX = clang
        CXX = g++
        TARGET = a.out
    endif
    ifeq ($(UNAME_S),Linux)
        CXX = g++
        TARGET = a.out
    endif
endif

CXXFLAGS = --stdlib=libc++ --std=c++0x -O3
SRCS = $(MAIN_CPP).cpp ibm.cpp sp_util.cpp
OBJS := $(SRCS:.cpp=.o)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS)

clean:
	rm -f $(TARGET) $(OBJS)