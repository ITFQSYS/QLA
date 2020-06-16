CXX := clang++
TARGET := test_main
INCLUDE := -I googletest/googletest/include/
LIBS := -L googletest/build/lib -lgtest -lgtest_main -lpthread
SOURCES = $(wildcard *.cpp) 


all: run_test

run_test : clean $(TARGET)
	./$(TARGET)

$(TARGET) : $(SOURCES) $(wildcard ../include/)
	$(CXX) $(SOURCES) -o $(TARGET) $(INCLUDE) $(LIBS)

.PHONY clean : 
	touch $(TARGET)
	rm  $(TARGET)