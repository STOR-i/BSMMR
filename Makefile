FLAGS=-std=c++0x -O3

BSMMR: BSMMR.o function_evaluation.o calculate_integral_birth.o modelfunctions.o calculate_accept_prob.o updates.o utilities.o cross_validation.o analysis.o
	g++ $^ -o $@

BSMMR.o: main.cpp
	g++ -c $(FLAGS) $< -o $@

%.o: %.cpp %.hpp
	g++ -c $(FLAGS) $< -o $@

clean:
	rm -f *.o test_function_evaluation BSMMR *-

