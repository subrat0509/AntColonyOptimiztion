/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package AntColonyOptimization_TSP;


import java.util.*;

public class AntColonySystemTSP implements AntOptimizationTSP{
    /* Distance Matrix */
    private int dist[][];
    /* Pheromone matrix */
    private double pheromone[][];
    /* matrix of heuristic values */
    private double choiceInfo[][];
    /* the vector of artificial ants */
    private SingleAnt ants[];

    /* number of ants */
    private int noAnts;
    /*  number of nodes */
    private int noNodes;
    /* alpha parameter */
    private double alfa;
    /* beta parameter  */
    private double beta;
    /* parameter ro = global evaporation*/
    private double ro;
    /*  parameter xi = local evaporation*/
    private double xi;
    /* initial amount of pheromones */
    private double tau0;
    /* the probability of choosing the shortest edge at a time  */
    private double q0;
    /* current iteration */
    private int iterations;
    /* if localsearch is applied  */
    private boolean doLocalSearch;
    /* if opt2 applies  */
    private boolean doOpt2;
    /*  if opt3 applies  */
    private boolean doOpt3;

    private OptimizationTSP opt;

    /*  the best tour found so far */
    private int[] bestSoFarTour;


    public AntColonySystemTSP(int noNodes, int noAnts){
        this.noNodes = noNodes;
        this.noAnts = noAnts;
        this.alfa = 1.0;
        this.beta = 3.0;
        this.ro = 0.1;
        this.xi = 0.1;
        this.iterations = 0;
        this.tau0 = 0;
        this.q0 = 0.9;

        dist = new int[noNodes][noNodes];
        pheromone = new double[noNodes][noNodes];
        choiceInfo = new double[noNodes][noNodes];

        ants = new SingleAnt[noAnts];
        for(int i = 0; i < noAnts; i++){
            ants[i] = new SingleAnt(noNodes);
        }
        doLocalSearch = false;
        bestSoFarTour = new int[noNodes+1];

    }
    public AntColonySystemTSP(int noNodes, int noAnts, double alfa,
                    double beta, double ro, double xi, double q0, boolean doOpt2, boolean doOpt3){
        this.noNodes = noNodes;
        this.noAnts = noAnts;
        this.alfa = alfa;
        this.beta = beta;
        this.ro = ro;
        this.xi = xi;
        this.iterations = 0;
        this.tau0 = 0;
        this.q0 = q0;

        dist = new int[noNodes][noNodes];
        pheromone = new double[noNodes][noNodes];
        choiceInfo = new double[noNodes][noNodes];

        ants = new SingleAnt[noAnts];
        for(int i = 0; i < noAnts; i++){
            ants[i] = new SingleAnt(noNodes);
        }
        //this.doLocalSearch = doLocalSearch;
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
    }
    public void setAlfa(double a){
        alfa = a;
    }
    public void setBeta(double b){
        beta = b;
    }
    public void setTau0(double tau){
        tau0 = tau;
    }
    public void setRo(double r){
        ro = r;
    }
    public void setXi(double x){
        xi = x;
    }
    public void setQ0(double q){
        q0 = q;
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
    public double getXi(){
        return xi;
    }
    public double getTau0(){
        return tau0;
    }
    public double getQ0(){
        return q0;
    }
    public int[] getBestTour(){
        int[] bestTour = new int[noNodes+1];
        int bestTourLength = Integer.MAX_VALUE;
        int bestIdx = -1;
        for(int i = 0; i < noAnts; i++){
            if(ants[i].getTourLength() < bestTourLength){
                bestTourLength = ants[i].getTourLength();
                bestIdx = i;
            }
        }
        for(int j = 0; j <=noNodes; j++)
            bestTour[j] = ants[bestIdx].getTour(j);
        return bestTour;
    }
    public int[] getBestSoFarTour(){
        return bestSoFarTour;
    }
    public void updateBestSoFarTour(){
        if(getBestTourLength() < computeTourLength(bestSoFarTour)){
            bestSoFarTour = getBestTour();
        }
    }
    public int getBestTourLength(){
        return computeTourLength(getBestTour());
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

        int c = ants[k].getTour(step-1); // the former city of the current ant
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
        /* The decision to be made if row <q0 */
        int randomDecision = j;

        /* the best decision */
        double maxHeuristic = -1;
        int maxHeuristicIdx = -1;
        for(j = 0; j < noNodes; j++){
            if(maxHeuristic < choiceInfo[c][j] && !(ants[k].getVisited(j))){
                maxHeuristic = choiceInfo[c][j];
                maxHeuristicIdx = j;
            }
        }

        if(Math.random() < q0){
            ants[k].setTour(step, maxHeuristicIdx);
            ants[k].setVisited(maxHeuristicIdx, true);
        }
        else{
            ants[k].setTour(step, randomDecision);
            ants[k].setVisited(randomDecision, true);
        }

        

    }
    public void constructSolutions(){
        /* erase memory ants */
        initAnts();

        int step = 0;
        int k;
        int r;

        Random rand = new Random();

        /* initial city assignment  */
        for(k = 0; k < noAnts; k++){
            r = Math.abs(rand.nextInt())%noNodes;

            ants[k].setTour(step,r);
            ants[k].setVisited(r,true);
        }
        /* effective construction of the solution  */
        while(step < noNodes-1){
            step++;
            for(k = 0; k < noAnts; k++){
                decisionRule(k,step);
                localPheromoneUpdate(k,step);
            }
        }
        /*  complete the tour  */
        for(k = 0; k < noAnts; k++){
            ants[k].setTour(noNodes,ants[k].getTour(0));
            localPheromoneUpdate(k,noNodes);
            ants[k].setTourLength(computeTourLength(ants[k].getTour()));
        }
        updateBestSoFarTour();
    }
    public void globalEvaporation(){
        /* the global evaporation takes place only on the arcs that belong to the best tour so far */
        for(int i = 0; i <noNodes; i++){
            int idx1 = bestSoFarTour[i];
            int idx2 = bestSoFarTour[i+1];
            //System.out.println("Pheromone before global evaporation depozit: "+pheromone[idx1][idx2]);
            pheromone[idx1][idx2]*=(1-ro);
            pheromone[idx2][idx1]*=(1-ro);
            //System.out.println("Pheromone after global evaporation depozit: "+pheromone[idx1][idx2]);
        }
        
    }
    public void depositPheromone(int k){
        /* storage only takes place on the best tour so far, so don't use this method */
        
    }
    public void globalPheromoneDeposit(){
        double delta = 1.0/((double) getBestTourLength());
        for(int i = 0; i <noNodes; i++){
            int idx1 = bestSoFarTour[i];
            int idx2 = bestSoFarTour[i+1];
            //System.out.println("Pheromone before global pheromone depozit: "+pheromone[idx1][idx2]);
            pheromone[idx1][idx2]+=ro*delta;
            pheromone[idx2][idx1]+=ro*delta;
            //System.out.println("Pheromone after global pheromone depozit: "+pheromone[idx1][idx2]);
        }
    }
    public void updatePheromones(){
        globalEvaporation();
        globalPheromoneDeposit();
        computeHeuristic();
    }
    public void localPheromoneUpdate(int ant, int step){
        int idx1 = ants[ant].getTour(step);
        int idx2 = ants[ant].getTour(step-1);
        //System.out.println("Pheromone before local evaporation:"+pheromone[idx1][idx2]);
        double currentValue = pheromone[idx1][idx2];
        pheromone[idx1][idx2] = (1-xi)*currentValue+xi*tau0;
        pheromone[idx2][idx1] = pheromone[idx1][idx2];
        //heuristic value update
        double niu = 0;
        if(dist[idx1][idx2] > 0)
            niu = 1.0/dist[idx1][idx2];
        else
            niu = 1.0/0.0001;
        //System.out.println("Pheromone after local evaporation:"+pheromone[idx1][idx2]);
        choiceInfo[idx1][idx2] = Math.pow(pheromone[idx1][idx2],alfa)*Math.pow(niu,beta);
        choiceInfo[idx1][idx2] = choiceInfo[idx2][idx1];
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
        return 1.0/(((double)greedyTour())*((double)noAnts));
    }
    public void opt2(){
        int i,j,k;
        int a1,a2,a3,b1,b2,b3,swap;
        /* for each ant */
        for(k = 0; k < noAnts; k++){
            //long len1 = computeTourLength(ants[k].getTour());
            for(i = 1; i < noNodes-1; i++){
                a1 = dist[ants[k].getTour(i-1)][ants[k].getTour(i)];
                a2 = dist[ants[k].getTour(i)][ants[k].getTour(i+1)];
                a3 = dist[ants[k].getTour(i+1)][ants[k].getTour(i+2)];

                b1 = dist[ants[k].getTour(i-1)][ants[k].getTour(i+1)];
                b2 = dist[ants[k].getTour(i+1)][ants[k].getTour(i)];
                b3 = dist[ants[k].getTour(i)][ants[k].getTour(i+2)];

                if(a1+a2+a3 > b1+b2+b3){
                    swap = ants[k].getTour(i);
                    ants[k].setTour(i, ants[k].getTour(i+1));
                    ants[k].setTour(i+1, swap);
                }

            }
            //long len2 = computeTourLength(ants[k].getTour());
            //if(len2 < len1)
            //    System.out.println("local search improvement: "+len1+" -> "+len2);
        }
    }
    public void localSearch(){
    /* Local search procedures */
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
    }
    public double[][] getPheromoneMatrix(){
        return pheromone;
    }


}
