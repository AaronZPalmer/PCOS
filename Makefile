
path="/Users/aaronpalmer/Documents/Examples/"

filename=PCOS.out

PCOS :
	c++ -Wall main.cpp PCOS.cpp PCOS.hpp Policy.cpp Policy.hpp PI.cpp PI.hpp -DFILEPATH=\"$(path)\"
	mv a.out $(filename)