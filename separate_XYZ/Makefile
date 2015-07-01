NAME := separate_XYZ  
EXTRA_SRC := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(EXTRA_SRC))

CC := g++
INCLUDE := -I/opt/local/include 
CPPFLAGS := -O3 -std=c++11 -lhdf5 -lhdf5_cpp #-fopenmp
LDFLAGS := -lhdf5 -lhdf5_cpp #-fopenmp 
LDLIBS := -L/opt/local/lib

$(NAME): $(OBJS)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o: %.c
	$(CC) $(CPPFLAGS) -c $<
clean:
	rm -f $(NAME) $(OBJS)

.PHONY: clean

