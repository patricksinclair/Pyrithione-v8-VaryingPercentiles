import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.ArrayList;
import java.util.Random;
import java.util.stream.IntStream;

class BioSystem {

    private Random rand = new Random();

    private double alpha, c_max; //steepness and max val of antimicrobial concn
    private double scale, sigma; //mic distb shape parameters
    private ArrayList<Microhabitat> microhabitats;
    private double time_elapsed, exit_time; //exit time is the time it took for the biofilm to reach the thickness limit, if it did
    private int immigration_index;

    private double deterioration_rate = 0.0516;
    private double immigration_rate = 0.8;
    private double tau = 0.01;
    private double delta_x = 5.;
    private int thickness_limit = 50; //this is how big the system can get before we exit. should reduce overall simulation duration
    private int n_detachments = 0, n_deaths = 0, n_replications = 0, n_immigrations = 0;


    private BioSystem(double alpha, double c_max, double scale, double sigma){

        this.alpha = alpha;
        this.c_max = c_max;
        this.scale = scale;
        this.sigma = sigma;
        this.microhabitats = new ArrayList<>();
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.immigration_index = 0;

        microhabitats.add(new Microhabitat(calc_C_i(0, this.c_max, this.alpha, delta_x), scale, sigma));
        microhabitats.get(0).setSurface();
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }

    private int getN_detachments(){return n_detachments;}
    private int getN_deaths(){return n_deaths;}
    private int getN_replications(){return n_replications;}
    private int getN_immigrations(){return n_immigrations;}

    private double getTimeElapsed(){return time_elapsed;}
    private double getExit_time(){return exit_time;}
    private int getBiofilmThickness(){return microhabitats.size();}

    private void setExit_time(double exit_time){this.exit_time = exit_time;}

    private int getTotalN(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN();
        }
        return runningTotal;
    }

    private int getBiofilmEdge(){
        int edgeIndex = 0;
        for(int i = 0; i < microhabitats.size(); i++){
            if(microhabitats.get(i).isBiofilm_region()) edgeIndex = i;
        }
        return edgeIndex;
    }




    private void immigrate(int mh_index, int n_immigrants){
        microhabitats.get(mh_index).addARandomBacterium_x_N(n_immigrants);
    }


    private void migrate(ArrayList<Microhabitat> microhabs, int mh_index, int bac_index){

        double migrating_bac = microhabs.get(mh_index).getPopulation().get(bac_index);
        microhabs.get(mh_index).removeABacterium(bac_index);

        if(microhabs.get(mh_index).isSurface()){
            microhabs.get(mh_index+1).addABacterium(migrating_bac);

        }else if(microhabs.get(mh_index).isImmigration_zone()){
            microhabs.get(mh_index-1).addABacterium(migrating_bac);

        }else{
            if(rand.nextBoolean()){
                microhabs.get(mh_index+1).addABacterium(migrating_bac);
            }else{
                microhabs.get(mh_index-1).addABacterium(migrating_bac);
            }
        }
    }



    private void updateBiofilmSize(){
        //once the edge microhabitat is sufficiently populated, this adds another microhabitat onto the system list
        //which is then used as the immigration zone

        if(microhabitats.get(immigration_index).atBiofilmThreshold()){

            microhabitats.get(immigration_index).setBiofilm_region();
            microhabitats.get(immigration_index).setImmigration_zone(false);

            int i = microhabitats.size();
            microhabitats.add(new Microhabitat(BioSystem.calc_C_i(i, c_max, alpha, delta_x), scale, sigma));
            immigration_index = i;
            microhabitats.get(immigration_index).setImmigration_zone(true);
        }

        //this stops sims going onn unnecessarily too long. if the biofilm reaches the thickness limit then we record the
        //time this happened at and move on
        if(getBiofilmThickness()==thickness_limit){
            exit_time = time_elapsed;
            time_elapsed = 9e9; //this way the time elapsed is now way above the duration value, so the simulation will stop
        }
    }


    private void performAction(){

        double tau_step = tau; //tau used for tau leaping time increment

        int system_size = microhabitats.size(); //this is all the microhabs in the system
        int[][] replication_allocations;
        int[][] death_allocations;
        int[][] migration_allocations;
        int[] detachment_allocations;
        int[] original_popsizes;
        int n_immigrants;

        whileloop:
        while(true){
            replication_allocations = new int[system_size][];
            death_allocations = new int[system_size][];
            migration_allocations = new int[system_size][];
            original_popsizes = new int[system_size];
            detachment_allocations = new int[microhabitats.get(immigration_index).getN()];

            for(int mh_index = 0; mh_index < system_size; mh_index++){
                int mh_pop = microhabitats.get(mh_index).getN();
                int[] n_replications = new int[mh_pop];
                int[] n_deaths = new int[mh_pop];
                int[] n_migrations = new int[mh_pop];

                for(int bac_index = 0; bac_index < mh_pop; bac_index++){

                    ////////// MIGRATIONS //////////////////////
                    n_migrations[bac_index] = new PoissonDistribution(microhabitats.get(mh_index).migrate_rate()*tau_step).sample();

                    if(n_migrations[bac_index] > 1){
                        tau_step /= 2.;
                        continue whileloop;
                    }
                    ////////////////////////////////////////////

                    ///////////// DETACHMENTS /////////////////////////
                    if(mh_index == immigration_index){
                        detachment_allocations[bac_index] = new PoissonDistribution(deterioration_rate*tau_step).sample();

                        if(detachment_allocations[bac_index] > 1){
                            tau_step /= 2.;
                            continue whileloop;
                        }
                        //if a bacteria is detaching then it can't migrate
                        if(detachment_allocations[bac_index] != 0){
                            n_migrations[bac_index] = 0;
                        }
                    }
                    ////////////////////////////////////////////////////////

                    ////////////////// REPLICATIONS AND DEATHS ///////////////////////////
                    double g_or_d_rate = microhabitats.get(mh_index).replicationOrDeathRate(bac_index);

                    if(g_or_d_rate == 0.){

                        n_replications[bac_index] = 0;
                        n_deaths[bac_index] = 0;

                    }else if(g_or_d_rate > 0){

                        n_replications[bac_index] = new PoissonDistribution(g_or_d_rate*tau_step).sample();
                        n_deaths[bac_index] = 0;

                    }else{
                        n_replications[bac_index] = 0;
                        n_deaths[bac_index] = new PoissonDistribution(Math.abs(g_or_d_rate)*tau_step).sample();

                        if(n_deaths[bac_index] > 1){
                            tau_step /= 2.;
                            continue whileloop;
                        }
                        //if a death is occurring, then that bacteria can't migrate or detach
                        if(n_deaths[bac_index] !=0) {
                            n_migrations[bac_index] = 0;
                            if(mh_index == immigration_index) detachment_allocations[bac_index] = 0;
                        }
                    }
                    /////////////////////////////////////////////////////////////////////////
                }

                replication_allocations[mh_index] = n_replications;
                death_allocations[mh_index] = n_deaths;
                migration_allocations[mh_index] = n_migrations;
                original_popsizes[mh_index] = microhabitats.get(mh_index).getN();
            }


            n_immigrants = new PoissonDistribution(immigration_rate*tau_step).sample();
            break;
        }


        for(int mh_index = 0; mh_index < system_size; mh_index++){
            for(int bac_index = original_popsizes[mh_index]-1; bac_index >= 0; bac_index--){

                if(death_allocations[mh_index][bac_index]!= 0) {
                    microhabitats.get(mh_index).removeABacterium(bac_index);
                    n_deaths++;
                }

                else{
                    microhabitats.get(mh_index).replicateABacterium_x_N(bac_index, replication_allocations[mh_index][bac_index]);
                    n_replications += replication_allocations[mh_index][bac_index];

                    if(system_size > 1){
                        if(migration_allocations[mh_index][bac_index] != 0) migrate(microhabitats, mh_index, bac_index);
                    }

                    if(mh_index == immigration_index){
                        if(detachment_allocations[bac_index] != 0) {
                            microhabitats.get(mh_index).removeABacterium(bac_index);
                            n_detachments++;
                        }
                    }
                }
            }
        }

        immigrate(immigration_index, n_immigrants);
        n_immigrations += n_immigrants;
        updateBiofilmSize();
        time_elapsed += tau_step;
    }



    private static int[] getThicknessAndEventCountersReachedAfterATime(double duration, int i, double scale, double sigma){
        int K = 120;
        double c_max = 10., alpha = 0.01;

        BioSystem bs = new BioSystem(alpha, c_max, scale, sigma);
        System.out.println("detach_rate: "+bs.deterioration_rate);
        int nUpdates = 20;
        double interval = duration/nUpdates;
        boolean alreadyRecorded = false;


        while(bs.time_elapsed <= duration){

            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int max_poss_pop = bs.getBiofilmThickness()*K;
                System.out.println("rep : "+i+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+bs.getTotalN()+"/"+max_poss_pop+"\tbf_edge: "+bs.getBiofilmEdge());
                alreadyRecorded = true;
            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;


            bs.performAction();
        }
        if((int)bs.exit_time == 0) bs.exit_time = duration;

        return new int[]{bs.getBiofilmEdge(), bs.getN_deaths(), bs.getN_detachments(), bs.getN_immigrations(), bs.getN_replications(), (int)bs.getExit_time()};
    }


    static void getBiofilmThicknessHistoInParallel(int nReps, double scale, double sigma){
        long startTime = System.currentTimeMillis();

        int nSections = 9; //number of sections the reps will be divided into, to avoid using loadsa resources
        int nRuns = nReps/nSections; //number of runs in each section

        double duration = 25.*7.*24.; //25 week duration

        int[][] index_and_counters_reached = new int[nReps][];

        String index_reached_filename = "pyrithione-t="+String.valueOf(duration)+"-parallel-event_counters_sigma="+String.format("%.5f", sigma);
        String[] headers = new String[]{"bf edge", "n_deaths", "n_detachments", "n_immigrations", "n_replications", "exit time"};

        for(int j = 0; j < nSections; j++){
            System.out.println("section: "+j);

            IntStream.range(j*nRuns, (j+1)*nRuns).parallel().forEach(i ->
                    index_and_counters_reached[i] = BioSystem.getThicknessAndEventCountersReachedAfterATime(duration, i, scale, sigma));
        }

        Toolbox.writeCountersToFile(index_reached_filename, headers, index_and_counters_reached);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }


    private static double calc_C_i(int i, double c_max, double alpha, double delta_x){
        return c_max*Math.exp(-alpha*i*delta_x);
    }

}
