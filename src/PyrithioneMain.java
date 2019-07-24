public class PyrithioneMain {
    public static void main(String[] args){

        double scale_99 = 2.71760274, sigma_99 = 0.56002833;
        double scale_95 = 1.9246899, sigma_95 = 1.00179994;
        double scale_90 = 1.01115312, sigma_90 = 1.51378016;
        String date = "-24-July-2019";
        String folderID99 = "-99_suscep"+date;
        String folderID95 = "-95_suscep"+date;
        String folderID90 = "-90_suscep"+date;

        int nReps = 180;
        System.out.println("geno distbs now with thickness limit and directory creation");

        BioSystem.getEventCountersAndRunPopulations(nReps, scale_99, sigma_99, folderID99);
    }
}
