/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package AntColonyOptimization_TSP;


import java.util.Random;

public class GeneticAntSystemTSP implements AntOptimizationTSP{
    /* distance matrix */
    private int dist[][];
    /* pheromone matrix  */
    private double pheromone[][];
    /* heuristic values ​​matrix  */
    private double choiceInfo[][];
    /* the vector of artificial ants */
    private SingleAnt ants[];

    /*number of ants  */
    private int noAnts;
    /*number of nodes */
    private int noNodes;
    /* alpha parameter */
    private double alfa;
    /* beta parameter  */
    private double beta;
    /* row parameter */
    private double ro;
    /* current iteration */
    private int iterations;
    /* if localsearch is applied */
    private boolean doLocalSearch;
    /* if opt2 applies */
    private boolean doOpt2;
    /* if opt3 applies */
    private boolean doOpt3;
    /* the best tour found so far */
    private int[] bestSoFarTour;
    private OptimizationTSP opt;

    private double tau0;

    /* min pheromone */
    private double minPheromone;
    /* max pheromone */
    private double maxPheromone;
    /* maximum number of iterations with stagnation  */
    private int maxStagnationIterations;
    /* current number of stagnant iterations  */
    private int stagnationIterations;
    /* number of iterations of the genetic algorithm */
    private int geneticIterations;

    /* applied genetic algorithm */
    private GeneticAlgorithmTSP ga;


    public GeneticAntSystemTSP(int noNodes, int noAnts, double alfa,
                    double beta, double ro,
                    boolean doOpt2, boolean doOpt3,
                    int maxStagnationIterations, int geneticIterations){
        this.noNodes = noNodes;
        this.noAnts = noAnts;
        this.alfa = alfa;
        this.beta = beta;
        this.ro = ro;
        this.iterations = 0;
        this.stagnationIterations = 0;
        this.maxStagnationIterations = maxStagnationIterations;
        this.geneticIterations = geneticIterations;

        dist = new int[noNodes][noNodes];
        pheromone = new double[noNodes][noNodes];
        choiceInfo = new double[noNodes][noNodes];

        ants = new SingleAnt[noAnts];
        for(int i = 0; i < noAnts; i++){
            ants[i] = new SingleAnt(noNodes);
        }

        ga = new GeneticAlgorithmTSP(noAnts, noNodes, 0.1, false, false);
        
        this.doOpt2 = doOpt2;
        this.doOpt3 = doOpt3;
        bestSoFarTour = new int[noNodes+1];
    }

    public SingleAnt getAnt(int k){
        return ants[k];
    }

    public double getPheromone(int i, int j){
        return pheromone[i][j];
    }
    public double getHeuristic(int i, int j){
        return choiceInfo[i][j];
    }
    public int getDistance(int i, int j){
        return dist[i][j];
    }
    public void setPheromone(int i, int j, double ph){
        pheromone[i][j] = ph;
    }
    public void setHeuristic(int i, int j, double h){
        choiceInfo[i][j] = h;
    }
    public void setDistance(int i, int j, int d){
        dist[i][j] = d;
        ga.setDistance(i, j, d);
    }
    public void setAlfa(double a){
        alfa = a;

    }
    public void setBeta(double b){
        beta = b;
    }
    public void setRo(double r){
        ro = r;
    }
    public double getAlfa(){
        return alfa;
    }
    public double getBeta(){
        return beta;
    }
    public double getRo(){
        return ro;
    }
    public int[] getBestTour(){
        int[] bestTour = new int[noNodes+1];
        int bestTourLength = Integer.MAX_VALUE;
        for(int i = 0; i < noAnts; i++){
            if(ants[i].getTourLength() < bestTourLength){
                bestTourLength = ants[i].getTourLength();
                for(int j = 0; j <=noNodes; j++)
                    bestTour[j] = ants[i].getTour(j);
            }
        }
        return bestTour;
    }
    public int getBestTourLength(){
        return computeTourLength(getBestTour());
    }
    public int[] getBestSoFarTour(){
        return bestSoFarTour;
    }
    public void updateBestSoFarTour(){
        if(getBestTourLength() < computeTourLength(bestSoFarTour)){
            bestSoFarTour = getBestTour();
            stagnationIterations = 0;
        }
        else
            stagnationIterations++;
    }
    public void setIteration(int iter){
        iterations = iter;
    }
    public int getIteration(){
        return iterations;
    }
    public int getNoAnts(){
        return noAnts;
    }
    public void setNoAnts(int ants){
        noAnts = ants;
    }
    public int getNoNodes(){
        return noNodes;
    }
    public void setNoNodes(int nodes){
        noNodes = nodes;
    }
    public void initData(){
        int i,j;
        for(i = 0; i < noNodes; i++)
            for(j = 0; j < noNodes; j++){
                dist[i][j] = 0;
                pheromone[i][j] = 0.0;
                choiceInfo[i][j] = 0.0;
            }
        ga.initPopulation();
        
    }
    public void initPheromones(){
        int i,j;
        tau0 = computePheromone0();

        for(i = 0; i < noNodes; i++)
            for(j = 0; j < noNodes; j++)
                pheromone[i][j] = tau0;

        for(i = 0; i < noNodes; i++)
            pheromone[i][i] = 0;

        opt = new OptimizationTSP(dist);

        maxPheromone = tau0;
        minPheromone = maxPheromone*(1-Math.pow(0.05, 1.0/((double)noNodes)))/(noNodes);

    }
    public void reinitPheromones(){
        int i,j;


        for(i = 0; i < noNodes; i++)
            for(j = 0; j < noNodes; j++)
                pheromone[i][j] = tau0;

        for(i = 0; i < noNodes; i++)
            pheromone[i][i] = 0;

        //opt = new OptimizationTSP(dist);

        maxPheromone = tau0;
        minPheromone = maxPheromone*(1-Math.pow(0.05, 1.0/((double)noNodes)))/(noNodes);

    }
    public void computeHeuristic(){
        double niu;
        int i,j;

        for(i = 0; i < noNodes; i++)
            for(j = 0; j < noNodes; j++){
                if(dist[i][j] > 0)
                    niu = 1.0/dist[i][j];
                else
                    niu = 1.0/0.0001;
            choiceInfo[i][j] = Math.pow(pheromone[i][j],alfa)*Math.pow(niu,beta);
        }
    }
    public void initAnts(){
        int i,j;
        for(i = 0; i < noAnts; i++){
            ants[i].setTourLength(0);

            for(j = 0; j < noNodes; j++)
                ants[i].setVisited(j, false);
            for(j = 0; j <= noNodes; j++)
                ants[i].setTour(j, 0);
        }
    }
    public void decisionRule(int k, int step){
        /* k = ant identifier */
        /* step = the current step in the construction of the solution */

        int c = ants[k].getTour(step-1); // the previous city of the current ant
        double sumProb = 0.0;

        double selectionProbability[] = new double[noNodes];

        int j;
        for(j = 0; j < noNodes; j++){
            if((ants[k].getVisited(j)) || (j == c))
                selectionProbability[j] = 0.0;
            else{
                selectionProbability[j] = choiceInfo[c][j];
                sumProb+=selectionProbability[j];
            }

        }
        double prob = Math.random()*sumProb;
        j = 0;
        double p = selectionProbability[j];
        while(p < prob){
            j++;
            p += selectionProbability[j];
        }
        ants[k].setTour(step, j);
        ants[k].setVisited(j, true);

    }
    public void constructSolutions(){
        /* erase memory ants  */
        initAnts();

        int step = 0;
        int i,k;
        int r;

        Random rand = new Random();

        /*  initial city assignment*/
        for(k = 0; k < noAnts; k++){
            r = Math.abs(rand.nextInt())%noNodes;

            ants[k].setTour(step,r);
            ants[k].setVisited(r,true);
        }
        /* effective construction of the solution */
        while(step < noNodes-1){
            step++;
            for(k = 0; k < noAnts; k++)
                decisionRule(k,step);
        }
        /* complete the tour */
        for(k = 0; k < noAnts; k++){
            ants[k].setTour(noNodes,ants[k].getTour(0));
            ants[k].setTourLength(computeTourLength(ants[k].getTour()));
        }
        updateBestSoFarTour();
        /* assigning the obtained solutions to the genetic algorithm*/
        for(k = 0; k < noAnts; k++){
            for(i = 0; i <= noNodes; i++){
                ga.getChromosome(k).setGene(i, ants[k].getTour(i));
            }
        }
        
        ga.updateBestSoFarTour();
        for(i = 0; i < geneticIterations; i++){
            ga.stepGeneration();
        }
        if(computeTourLength(bestSoFarTour) > computeTourLength(ga.getBestSoFarTour())){
            System.out.println("iddqd");
            int []newTour = ga.getBestSoFarTour();
            for(i = 0; i < noNodes+1; i++)
                bestSoFarTour[i] = newTour[i];
        }
        
        updateBestSoFarTour();
        
    }
    public void globalEvaporation(){
        int i,j;
        for(i = 0; i < noNodes; i++)
            for(j = 0; j < noNodes; j++){
                pheromone[i][j] = (1-ro)*pheromone[i][j];
            }
    }
    public void depositPheromone(int k){
        // Aici nu fac nimic

    }
    public void globalPheromoneDeposit(){
        double delta;
        if(iterations%2 == 0)
            delta = 1.0/((double) getBestTourLength());
        else
            delta = 1.0/((double) computeTourLength(getBestTour()));
        for(int i = 0; i <noNodes; i++){
            int idx1 = bestSoFarTour[i];
            int idx2 = bestSoFarTour[i+1];
            //System.out.println("Pheromone before global pheromone depozit: "+pheromone[idx1][idx2]);
            if(delta < minPheromone)
                delta = minPheromone;
            if(delta > maxPheromone)
                delta = minPheromone;
            pheromone[idx1][idx2]+=delta;
            pheromone[idx2][idx1]+=delta;
            //System.out.println("Pheromone after global pheromone depozit: "+pheromone[idx1][idx2]);
        }
    }
    public void updatePheromones(){
        globalEvaporation();
        //for(int k = 0; k < noAnts; k++)
        //    depositPheromone(k);
        globalPheromoneDeposit();

        computeHeuristic();
        maxPheromone = 1/(ro*((double)computeTourLength(getBestSoFarTour())));
        minPheromone = maxPheromone*(1-Math.pow(0.05, 1.0/((double)noNodes)))/(noNodes);
        if(stagnationIterations > maxStagnationIterations){
            System.out.println("Stagnation Achieved");
            reinitPheromones();
            stagnationIterations = 0;
        }
    }
    private int greedyTour(){
        boolean visited[] = new boolean[noNodes];
        int tour[] = new int[noNodes+1];
        int length;
        int min, node;
        int i,j;

        for(i = 0; i < noNodes; i++)
            visited[i] = false;

        tour[0] = 0;
        bestSoFarTour[0] = 0;
        visited[0] = true;

        for(i = 1; i < noNodes; i++){
            min = Integer.MAX_VALUE;
            node = -1;
            for(j = 0; j < noNodes; j++){
                if((!visited[j])&&(j!=tour[i-1])){
                    if(min > dist[tour[i-1]][j]){
                        min = dist[tour[i-1]][j];
                        node = j;
                    }
                }
            }
            tour[i] = node;
            bestSoFarTour[i] = node;
            visited[node] = true;
        }
        tour[noNodes] = tour[0];
        bestSoFarTour[noNodes] = bestSoFarTour[0];
        return computeTourLength(tour);

    }
    public int computeTourLength(int tour[]){
        int len = 0;
        for(int i = 0; i < noNodes; i++){
            len+=dist[tour[i]][tour[i+1]];
        }
        return len;
    }
    private double computePheromone0(){
        return 1.0/(((double)greedyTour())*ro);
    }
    
    public void localSearch(){
    
        if(doOpt2){
            for(int i = 0; i < noAnts; i++){
                opt.opt2(ants[i].getTour());
            }
        }
        if(doOpt3){
            for(int i = 0; i < noAnts; i++){
                opt.opt3(ants[i].getTour());
            }
        }


        updateBestSoFarTour();
        if(stagnationIterations > 0)
            stagnationIterations--;
    }
    @Override
    public String toString(){
        String str = "";
        for(int i = 0; i < noAnts; i++){
            str = str+"{"+i+"}: ";
            for(int j = 0; j < noNodes+1; j++){
                str+=""+ants[i].getTour(j)+", ";
            }
            str+="[ "+computeTourLength(ants[i].getTour())+" ]";
            str+="\n";
        }
        return str;
    }
    public double[][] getPheromoneMatrix(){
        return pheromone;
    }

}

