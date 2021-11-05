
import java.lang.Math;

public class Bayes {
    double[] data = new double[0];
    public Bayes(double[] data){
        this.data = data;
    }
    // pMM 
    public double get5mC5mC(){
        double pi = 1;
        for(int i = 0; i<data.length; i++){
            double p = data[i]/2.0 + data[i]/2.0;
            pi *= p;
        }
        if (pi == 0){
            pi = 1.0E-12;
        }
        return -10* (Math.log10(0.04040897) + Math.log10(pi));
    }
    // pMC 
    public double get5mCC(){
        double pi = 1;
        for(int i = 0; i<data.length; i++){
            double p = data[i]/2.0 + (1.0-data[i])/2.0;
            pi *= p;
        }
        if (pi == 0){
            pi = 1.0E-12;
        }
        return -10* (Math.log10(0.01794966) + Math.log10(pi));
    }
    // pCC 
    public double getCC(){
        double pi = 1;
        for(int i = 0; i<data.length; i++){
            double p = (1.0-data[i])/2.0 + (1.0-data[i])/2.0;
            pi *= p;
        }
        if (pi == 0){
            pi = 1.0E-12;
        }
        return -10* (Math.log10(0.9416414) + Math.log10(pi));
    }
}
