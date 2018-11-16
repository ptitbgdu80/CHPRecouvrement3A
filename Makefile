# Compilateur utilisé
CC=mpic++

run : main.cc fannyestlaplusjolie.cpp fannyestlaplusjolie.h DataFile.cpp DataFile.h
	$(CC) -std=c++11 -I Eigen/Eigen main.cc fannyestlaplusjolie.cpp DataFile.cpp -o run -Wall

#si on a des trucs a tester :
test : test.cc
	$(CC) test.cc  $(CXX_FLAGS) -o run_test

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ run
