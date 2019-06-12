import org.apache.commons.math3.distribution.LogNormalDistribution;

import java.util.ArrayList;


class Microhabitat {

    private LogNormalDistribution MIC_distribution;

    private double c; //concn of antimicrobial
    private ArrayList<Double> population; //list of MICs of bacteria in microhab

    private int K = 120; //karryking kapacity
    private boolean surface = false, biofilm_region, immigration_zone = false;

    Microhabitat(double c, double scale, double sigma){
        this.c = c;
        double mu = Math.log(scale);
        this.population = new ArrayList<>(K);

        this.MIC_distribution = new LogNormalDistribution(mu, sigma);
        biofilm_region = false;
    }



    int getN(){return population.size();}
    boolean isSurface(){return surface;}
    boolean isBiofilm_region(){return biofilm_region;}
    boolean isImmigration_zone(){return immigration_zone;}
    ArrayList<Double> getPopulation(){return population;}

    void setSurface(){this.surface = true;}
    void setBiofilm_region(){this.biofilm_region = true;}
    void setImmigration_zone(boolean immigration_zone){this.immigration_zone = immigration_zone;}


    private double fractionFull(){
        return getN()/(double)K;
    }

    boolean atBiofilmThreshold(){
        double threshold_density = 0.8;
        return fractionFull() >= threshold_density;}

    double migrate_rate(){
        //returns 0.5*b for the microhabitat next to the ship hull, to account for the inability to migrate into the hull
        //also for the microhabitat that's the biofilm edge
        double b = 0.2;
        return (surface || immigration_zone) ? 0.5*b : b;
    }

    private double beta(int index){
        return population.get(index);
    }

    private double phi_c(int index){
        double cB = c/beta(index);
        return 1. - (6.*cB*cB)/(5. + cB*cB);
    }

    double replicationOrDeathRate(int index){
        double phi_c_scaled = 0.083*(phi_c(index));
        return (phi_c(index) > 0.) ? phi_c_scaled*(1. - getN()/(double)K) : phi_c_scaled;
    }


    void addARandomBacterium_x_N(int n_bacteria){
        for(int i = 0; i < n_bacteria; i++){
            population.add(MIC_distribution.sample());
        }
    }

    void replicateABacterium_x_N(int index, int nReps){
        for(int i = 0; i < nReps; i++){
            population.add(population.get(index));
        }
    }

    void addABacterium(double MIC){population.add(MIC);}

    void removeABacterium(int index){population.remove(index);}



}
