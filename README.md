# code_multi_period
périodisation, avec ripmap et rotations...
Pour avoir methode naive : taper 0
pour avoir la methode avec szeliski non périodisé : taper 1
pour avoir la décomposition avec szeliski et le ripmap : taper 2
pour avoir la décomposition avec yaroslavsky et le ripmap : taper 3
pour avoir une image blanche : taper4 -_-
pout avoir un ripmap : taper 5


Pour la compilation :
c99 -O2 -DNDEBUG viho.c iio.c ftr.c -lX11 -lpng -ljpeg -ltiff -lfftw3 -o viho -lm  
