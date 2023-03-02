target 	:= main 
objs 	:= ./build/main.o ./build/functions.o ./build/matrices.o ./build/plotting.o ./build/integration.o ./build/RungeKutta.o ./build/fourier.o 
CC 		:= clang++
CFLAGS 	:= -O3 -fsanitize=address
# ASAN := verbosity=1:malloc_context_size=20 

all: $(target)
run: $(target)
	./$(target)
deps := $(patsubst %.o, %.d, $(objs))
# -include $(deps)
DEPFLAGS = -MMD -MF $(@:.o=.d)

main: $(objs)
	$(CC) $(CFLAGS) -o $@ $^

build/%.o: ./MagneticDynamics/%.cpp	
	$(CC) $(CFLAGS) -c  $<  -o $@ $(DEPFLAGS)

clean:
	rm -f $(target) $(objs) $(deps)