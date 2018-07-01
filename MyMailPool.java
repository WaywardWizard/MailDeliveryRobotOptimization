package strategies;

import automail.MailItem;
import automail.StorageTube;
import exceptions.TubeFullException;
import automail.Building;
import automail.Clock;
import automail.PriorityMailItem;

import java.util.ArrayList;
import java.util.List;
import java.util.Comparator;
import java.util.stream.Collectors;

public class MyMailPool implements IMailPool {
    
    /* Priority level of standard mail items as per the specification */
    private static int DEFAULT_PRIORITY = 0;
    private static double TIME_MEASURE_POWER_FACTOR = 1.1;
    private static int WEAK_ROBOT_ITEM_WEIGHT_LIMIT = 2000;
    private static int BOT_COUNT=2;
    private static int MAX_CARRY=4;
    
    /* A larger pool size gives mail loading more ability to minimize distances
     * travelled however increases time.
     */
    private static int CANDIDATE_POOL_SIZE=BOT_COUNT*MAX_CARRY;  // Defined as size total robot fleet carry capacity
    private static int CANDIDATE_MAIL_SCOPE=10; 
    private static int N_MAIL_CLUSTER=1;
 
    /***
     * The mail pool uses these to estimate when a bot will return to refill.
     * This information is used by the mail pool to optimize loads across bots.
     * 
     * A value of -1 indicates return time unknown. An integer value represents
     * the time estimated to return.
     */
    private int weakBotReturn = -1;
    private int strongBotReturn = -1;
    
    private ArrayList<MailItem> candidateMailPool = new ArrayList<MailItem>();
    private ArrayList<ArrayList<MailItem>> candidateMailCluster 
        = new ArrayList<ArrayList<MailItem>>();
    
    public MyMailPool() {
    }
    
    /***
     * Mail age comparator for priority queue sorting by age
     * @author ben
     */
    private class MailAgeComparator implements Comparator<MailItem> {
        public int compare(MailItem a, MailItem b) {
            return(a.getArrivalTime()-b.getArrivalTime());
        }
    }
    
    /***
     * Mail destination floor comparator. Used for finding mail item with maximum
     * destination floor among a group of mail items. 
     * @author ben
     *
     */
    private class MailDestinationFloorComparator implements Comparator<MailItem> {
        public int compare(MailItem a, MailItem b) {
            return(a.getDestFloor()-b.getDestFloor());
        }
    }
    
    /***
     * Comparator for low to high cost ordering of paths 
     * @author ben
     *
     */
    private class PathCostComparator implements Comparator<List<MailItem>> {
        public int compare(List<MailItem> a, List<MailItem> b) {
            Path aPath = new Path(a);
            Path bPath = new Path(b);
            double delta = aPath.getDeliveryRunCost()-bPath.getDeliveryRunCost(); 
            return((int)(Math.signum(delta)));
        }
    }

    /***
     * Utility object for path analysis
     * 
     * @author ben
     *
     */
    private class Path {
        
        /* Time costs are as follows; */
        private int DELIVERY_TCOST=1;
        private int PACK_FILL_TCOST=1;
        private int FLOOR_CHANGE_COST=1;
        private List<MailItem> path;
        
        /* Memoize paths costs for efficiency. Each element represents the 
         * time to make the delivery since departing the mailroom.  */
        private int[] legTime;
        private int[] cumulativeLegTime;
        private double[] deliveryLegCost;
        private double deliveryRunCost = 0;

        /***
         * Dummy constructor for empty paths 
         */
        public Path() {}
        
        public Path(List<MailItem> p) {
            path=p;

			legTime = new int[path.size()];
			cumulativeLegTime = new int [path.size()];
			deliveryLegCost = new double[path.size()];

            /* Preprocess leg times */
            int legTimeSummand=0;
            int legTimeLength;
            double legCost;
            for(int ix=0;ix<path.size();ix++) {
                legTimeLength = pathLeg(ix);
                legTimeSummand+=legTimeLength;
                legTime[ix]=legTimeLength;
                cumulativeLegTime[ix]=legTimeSummand;
                legCost = costLeg(ix);
                deliveryLegCost[ix]=legCost;
                deliveryRunCost+=legCost;
            }
        }
        
        /***
         * Return the time elapsed in a given leg of the path. The bot journey 
         * starts at the mail room. The first leg includes one unit for refilling
         * the tube. A delivery takes one unit of time. 
         * 
         * A leg is typically then;
         *      floorDistance + 1 
         * @param leg
         */
        private int pathLeg(int leg) {
            int currentFloor = 
                    (leg==0) ? Building.MAILROOM_LOCATION : path.get(leg-1).getDestFloor();
            int fillPackTime = (leg==0)?PACK_FILL_TCOST:0;
            int floorTraverseTime = FLOOR_CHANGE_COST * Math.abs(path.get(leg).getDestFloor()-currentFloor);
            int deliveryTime = DELIVERY_TCOST;
            return(fillPackTime + floorTraverseTime + deliveryTime);
        }
        
        /*** 
         * Returns true if this is an empty path (ie no mail to be delivered)
         * @return
         */
        public boolean pathEmpty() {
        	return(path==null);
        }
                
 
        private double getDeliveryRunCost() { return(deliveryRunCost); }

        /**
         * This function calculates the cost of delivering a mail item in a given
         * leg of this path. The cost is calculated by the cost function given
         * in the assignment spec; 
         * 
         *      Td - Mail delivery time
         *      Ta - Arrival time
         *      l  - Time power factor
         *      Mp - Mail priority
         *      
         *      Cost = (Td - Ta)^l x (1+sqrt(Mp))
         *      
         * Cost a leg of the robot's journey.
         * @param leg - the leg of the robot's journey
         * @return
         */
        private double costLeg(int leg) {
            int timeDelta = ((Clock.Time() + cumulativeLegTime[leg])-path.get(leg).getArrivalTime());
            double timeFactor = Math.pow(timeDelta, TIME_MEASURE_POWER_FACTOR);
            double priorityFactor = 1 + Math.sqrt((double)resolvePriority(path.get(leg)));
            return(timeFactor * priorityFactor);
        }

        
        public List<MailItem> getPathList() {return(path);}

        public int getReturnTime() {
        	/* Return time is immediate if there are no deliverys in path */
        	if (pathEmpty()) {return(Clock.Time());}

        	int homeTime = Math.abs(path.get(-1+path.size()).getDestFloor() - Building.MAILROOM_LOCATION);
            int currentTime = Clock.Time();
        	int deliveryTime = cumulativeLegTime[-1+cumulativeLegTime.length];
        	return(homeTime + currentTime + deliveryTime);
        }

        /***
         * Resolves the priority of a given MailItem. If the MailItem is a 
         * PriorityMailItem subclass then use its assigned priority value. Otherwise
         * use the default mail item priority.
         * @param a - Mail item for which to resolve priority
         * @return
         */
        private int resolvePriority(MailItem a) {
            if (a instanceof PriorityMailItem) {
                return(((PriorityMailItem) a).getPriorityLevel());
            } else {
                return DEFAULT_PRIORITY;
            }
        }
        
    }
    
    private MailAgeComparator ageComparator = new MailAgeComparator();
    private MailDestinationFloorComparator destComparator = new MailDestinationFloorComparator();
    private PathCostComparator costComparator = new PathCostComparator();
    private ArrayList<MailItem> priorityMailPool = new ArrayList<MailItem>();
    private ArrayList<MailItem> mailPool = new ArrayList<MailItem> ();

    @Override
    /**
     * Adds a mail item into the mail pool. 
     */
    public void addToPool(MailItem mailItem) {
        if (mailItem instanceof PriorityMailItem) {
            priorityMailPool.add(mailItem);
            priorityMailPool.sort(ageComparator);
        } else {
            mailPool.add(mailItem);
            mailPool.sort(ageComparator);
        }
    }

    @Override
    /**
     * Fill a storage tube with mail from the pool.
     */
    public void fillStorageTube(StorageTube tube, boolean strong) {

        clearCandidateMailPool();
        
        /* Travel distance optimized clusters */
        getCandidatePriorityMail();
        getCandidateMail();
        
        Path bestPlan = optimizeDeliveryPlan(strong);
        
        /* Short circuit if no mail to deliver */
        if (bestPlan.getPathList()==null) {return;}
     	

		/* Update bot return time */
		if (strong) {
			strongBotReturn=bestPlan.getReturnTime();
		} else {
			weakBotReturn = bestPlan.getReturnTime();
		}
   
		/* Tube is a stack - fill it backwards for correct removal order*/
        MailItem currentMail;
		for(int ix=(bestPlan.getPathList().size()-1);ix>=0;ix--) {
			try{
				
				/* Load mail into robot backpack and remove from mail pool */
				currentMail=bestPlan.getPathList().get(ix);
				tube.addItem(currentMail);
				priorityMailPool.remove(currentMail);
				mailPool.remove(currentMail);

			} catch (TubeFullException e) {
				break;
			}
		}
	}
    
    /***
     * For each possible set of delivery items the refilling bot could carry, called a
     * 'carry set' (taken from the candidate delivery pool), find the set of 
     * ordered permutations for that carry set (where each permutation
     * comprises a delivery path) and then find the permutation with the least
     * cost. 
     * 
     * The complimentary 'carry set' is the set which the next robot to return
     * will carry. It is assumed that the next robot to return is not the current
     * robot. In the same manner as above, score the complimentary carry set permutations
     * for the next bot, based on it's predicted return time, based on the estimated
     * time cost for given actions. 
     * 
     * Given the two costs of the optimal carry set permutations, sum them to 
     * get the 'plan cost'. 
     * 
     * Search across all possible carry sets for the lowest possible plan cost. 
     * 
     * Return the lowest possible plan cost. This includes;
     *      -> Carry Set for refilling bot, ordered in the form of it's optimal permutation
     *      
     * Note that the candidate mail pool should be initialized before calling this
     * function.
     * 
     * @param strong - indicates wether search for strong or weak robot
     * @return - Optimal path for strong or weak robot (depending on arg strong)
     * Return null if no mail ready for delivery.
     */
    public Path optimizeDeliveryPlan(boolean strong) {
    	if(candidateMailPool.size()==0) {return(new Path());}
        ArrayList<List<MailItem>> weakCarrySet = getListWeakCarrySet();
        ArrayList<List<MailItem>> strongCarrySet = getListComplimentaryCarrySet(weakCarrySet);
        
        Path optimalWeakPath=null;
        Path localOptimalWeakPath;
        Path optimalStrongPath=null;
        Path localOptimalStrongPath;
        
        double bestPlanCost = -1;
        double planCost=0;
        
        /* Find optimal plan cost for each carry set */
        for(int ix=0;ix<weakCarrySet.size();ix++){
            
        	planCost=0;
        	
            /* Find optimal carry set permutes*/
            localOptimalWeakPath = getOptimalPath(weakCarrySet.get(ix));
            localOptimalStrongPath = getOptimalPath(strongCarrySet.get(ix));

            /* Cost plan */
            planCost += localOptimalWeakPath.getDeliveryRunCost();
            planCost += localOptimalStrongPath.getDeliveryRunCost();
            
            /* Find minimum plan cost */
            if (bestPlanCost<0 || planCost<bestPlanCost) {
                optimalWeakPath = localOptimalWeakPath;
                optimalStrongPath = localOptimalStrongPath;
                bestPlanCost=planCost;
            }
        }
        
        Path plan=strong?optimalStrongPath:optimalWeakPath;
        return(plan);
    }


    /***
     * Return a list of possible carry sets assembleable from candidate pool.
     * 
     * All possible sets from the candidate pool for the weak bot will not contain 
     * any heavy items. 
     */
    private ArrayList<List<MailItem>> getListWeakCarrySet() {
        /* Set the minimum number of items the weak bot should carry to the most
         * the strong bot could leave over. If the strong bot might leave over
         * nothing, then the weak bot won't be required to carry items.*/
    	int minimumCarrySize = Math.max(0, (MAX_CARRY - getCandidateMailSpaces()));

    	/* Return all possible carry sets */
        return(setPermutator(minimumCarrySize, MAX_CARRY, candidateMailPool.stream()
					.filter(m->m.getWeight()<WEAK_ROBOT_ITEM_WEIGHT_LIMIT)
					.collect(Collectors.toList())));
    }
    

    /***
     * Given a list with each element, and List of mail items representing a 
     * 'carry set', find the compliment for each carry set and return it in a 
     * structure the same as the argument. 
     * @param l
     * @return
     */
    private ArrayList<List<MailItem>> 
    getListComplimentaryCarrySet(ArrayList<List<MailItem>> l) {

        ArrayList<List<MailItem>> complimentarySet 
            = new ArrayList<List<MailItem>>();
        
        List<MailItem> compliment;
        
        for(List<MailItem> set:l) {

            /* The compliment will always have a maximum of four items since the
             * weak carry set has a minimum carry set size that ensures this */
			compliment = (ArrayList<MailItem>)(candidateMailPool.clone());
			compliment.removeAll(set);

            /* Get compliment */
            complimentarySet.add(new ArrayList<MailItem>(compliment));
        }
        
        return(complimentarySet);
    }
    
    
    /***
     * Optimize cluster delivery for a carry set by brute force searching the 
     * solution space for the optimal solution. 
     * 
     * Return null if no delivery. (set is empty) Else return the optimal path
     */
    private Path getOptimalPath (List<MailItem> set) {
    	if(set.size()==0) {return(new Path());}
        ArrayList<List<MailItem>> permutations = combinatronicPermutator(set);
        return(new Path((permutations.stream().min(costComparator)).get()));
    }
    
    
    /***
     * Returns an array list of all possible permutations of a mail item ordering. 
     * Where n is the number of mail items in the generating ArrayList<MailItem>
     * This takes n! time. However because n = 4, this hardly amounts to anything.
     * 
     * @param c - Generating array list of mail items.
     * @return - n! entries, each representing a permutation of the argument
     */
    private ArrayList<List<MailItem>> combinatronicPermutator(List<MailItem> c) {

        ArrayList<List<MailItem>> permutation = new ArrayList<List<MailItem>>();

        /* Recursive base cases */
        if (c.size()<1) {return(permutation);}
        if (c.size()==1) {
            permutation.add(c);
            return(permutation);
        }

        /* Get permutations via recursion */
        for(MailItem mx:c) {
            
            /* Get a list of mail items to recurse on */
            List<MailItem> recurseList = new ArrayList<MailItem>(c);
            recurseList.remove(mx);

            /* Get permutations for n-1 list (aka recurseList) */
            ArrayList<List<MailItem>> priorPermutation = combinatronicPermutator(recurseList);

            /* Prefix for nth level permutations */
            ArrayList<MailItem> prefix = new ArrayList<MailItem>();
            prefix.add(mx);

            /* Combine prefix with recursion permutations */
            for(int ix=0;ix<(factorial(c.size()-1));ix++) {
                List<MailItem> permute = new ArrayList<MailItem>(prefix);
                permute.addAll(priorPermutation.get(ix));
                permutation.add(permute);
            }
        }

        return permutation;
    }
 

    /***
     * Returns all unique sets (unordered by definition) of sizes between and including 
     * 'minSize' and 'maxSize' assembled from source. 
     * 
     * If the source has less unique elements in it than specified by 'setSize' then
     * return the source. 
     * 
     * @param minSize
     * @param maxSize
     * @param source
     * @return - List of all subsets. Empty if size source < minSize
     */
    private ArrayList<List<MailItem>> setPermutator(int minSize, int maxSize, List<MailItem>source) {
    	
        /* The source should be a set (unique elements only)  */
        source = source.stream().distinct().collect(Collectors.toList());
        ArrayList<List<MailItem>> setPermute = new ArrayList<List<MailItem>>();

    	/* If the source size is less than the min size, no permutation exists */
    	if (source.size()<minSize) {
    		return(setPermute);
    	}

    	/* Find all sets of in the set size range that can be made from the source set */
    	for(int ix=minSize;ix<=maxSize;ix++) {
    		if (source.size()>=ix) {
				setPermutatorRec(source, ix, 0, new ArrayList<MailItem>(), setPermute);
    		}
    	}
        
        return(setPermute);
    }
    

    /** 
     * Recursive algorithm for finding all possible subsets of superset
     * Recursively checks all possibilities
     * 
     * @param superset
     * @param size -subset size
     * @param depthx - depth through source set 
     * @param subset- subset being constructed 
     * @param setList - list to add subsets found
     */
    private void setPermutatorRec(List<MailItem> superset, 
    	int size,
    	int depthx,
    	List<MailItem> subset,
    	ArrayList<List<MailItem>> setList) {

    	/* Base case - made up set of required size */
    	if (subset.size()==size) {
    		setList.add(new ArrayList<>(subset));
    		return;
		}

    	/* Base cases - cant make up set of required size from assumptions */
    	if (depthx>=superset.size()) {return;}

    	/* Assume in set and recurse */
    	subset.add(superset.get(depthx));
    	setPermutatorRec(superset, size, (depthx+1), subset, setList);

    	/* Assume not in set and recurse */
    	subset.remove(superset.get(depthx));
    	setPermutatorRec(superset, size, (depthx+1), subset, setList);
    	return;
     }
    

    private int factorial(int n) {
        if (n<=1) {return(1);}
        return(n*factorial(n-1));
    }


    /***
     * Empty the candidate mail pool
     */
    private void clearCandidateMailPool() {
        candidateMailPool.clear();
    }


    /***
     * Empty the candidate clusters
     */
    private void clearCandidateCluster() {
        for(ArrayList<MailItem> c:candidateMailCluster) {
            c.clear();
        }
    }
    

    /***
     * Clusters mail by delivery destination. Divide candidate mail into n
     * groups. Minimize the distance between the group members.
     * 
     * Clusters seeded from the farthest two mail items. This is an alternative 
     * to brute force searching the entire solution space & can be used for mailbots
     * with larger carry capacity. Its not used in this solution though.
     */
    private void clusterCandidateMail() {
        clearCandidateCluster();

        /* Seed clustering algorithm */
        MailItem maxSeed = candidateMailPool.stream().max(destComparator).get();
        MailItem minSeed = candidateMailPool.stream().min(destComparator).get();
        
        ArrayList<MailItem> seed = new ArrayList<MailItem>();
        seed.add(maxSeed);
        seed.add(minSeed);

        /* Initialize clusters with seeds */
        candidateMailPool.remove(maxSeed);
        candidateMailPool.remove(minSeed);
        candidateMailCluster.get(0).add(maxSeed);
        candidateMailCluster.get(1).add(minSeed);

        int r=0;

        /* Populate mail clusters from candidate mail pool */
        while (candidateMailPool.size()>0) {

            /* Grow cluster from seed */
            for (int cx=0;cx<N_MAIL_CLUSTER;cx++) {
            	
            	/* These need to be final for use in lambda fn.*/
            	final int destinationFloor = seed.get(cx).getDestFloor();
            	final int rad = r;

                List<MailItem> convert = candidateMailPool.stream()
                    .filter(m -> m.getDestFloor() < (destinationFloor+rad+1))
                    .filter(m -> m.getDestFloor() > (destinationFloor-rad-1))
                    .collect(Collectors.toList());

                candidateMailCluster.get(cx).addAll(convert);
                candidateMailPool.removeAll(convert);

            }
            
            /* Expand radius */
            r++;
        }
        
    }
    
    
    /***
     * Loads mail priority mail pool into candidate mail pool.  If the robot
     * is weak, dont give it heavy mail. If the robot is strong, prioritize all the
     * heavy mail available first and then the lightweight. Take as much priority
     * mail as can be carried.
     * 
     * Use any age ail items when collecting priority
     * Heavy mail is prioritized for strong bot. 
     * Elder mail is prioritized.
     * 
     * @param mailPool - pool from which to populate candidateMailPool
     * @param strong - if the robot whom will be taking the mail is strong
     */
    private void getCandidatePriorityMail() {

        /* Get heaviest of eldest mail. Dont get more than strong bot carry capacity,
         * since the remaining mail in the candidate pool will be optimized 
         * against the assumption that the weak bot will be carrying. */
        candidateMailPool.addAll(
             priorityMailPool.stream()
                 .filter(m -> m.getWeight()>=WEAK_ROBOT_ITEM_WEIGHT_LIMIT)
                 .limit(getCandidateHeavyMailSpaces())
                 .collect(Collectors.toList()));
        
        /* Get any remaining priority mail */
        candidateMailPool.addAll(
                priorityMailPool.stream()
                    // We already took the heavy mail
                     .filter(m->m.getWeight()<WEAK_ROBOT_ITEM_WEIGHT_LIMIT)
                     .limit(getCandidateMailSpaces()) 
                     .collect(Collectors.toList()));
    }
    
    
    private void getCandidateMail() {
        if(getCandidateMailSpaces()==0) {return;}

        /* Get heaviest of eldest mail. Dont get more than strong bot carry capacity,
         * since the remaining mail in the candidate pool will be optimized 
         * against the assumption that the weak bot will be carrying it. */
        candidateMailPool.addAll(
             mailPool.stream()
                 .limit(CANDIDATE_MAIL_SCOPE) // Dont worry about the recent mail. 
                 .filter(m -> m.getWeight()>=WEAK_ROBOT_ITEM_WEIGHT_LIMIT)
                 .limit(getCandidateHeavyMailSpaces())
                 .collect(Collectors.toList()));
    
        /* Get any remaining mail */
        candidateMailPool.addAll(
                mailPool.stream()
                     // We already took the heavy mail
                     .filter(m->m.getWeight()<WEAK_ROBOT_ITEM_WEIGHT_LIMIT)
                     .limit(getCandidateMailSpaces()) 
                     .collect(Collectors.toList()));
    }

    /***
     * Returns how much many candidate mail slots remain in the candidate mail 
     * pool. 
     * @return
     */
    private int getCandidateMailSpaces() {
        return(CANDIDATE_POOL_SIZE-candidateMailPool.size());
    }

    
    /***
     * Returns how many spaces there are for heavy mail in the candidate mail
     * pool. (this is determined by how many items the strong bot can carry.
     * @return - return 0 iif the candidate mail is full or if the heavy item
     * carry capacity is taken up.
     */
    private int getCandidateHeavyMailSpaces() {
        if(getCandidateMailSpaces()==0) {return(0);}
        int heavyMailCount = (int)(candidateMailPool.stream().filter(w->w.getWeight()>=WEAK_ROBOT_ITEM_WEIGHT_LIMIT).count());
        int space = MAX_CARRY-heavyMailCount;
        space = (space < 0)?0:space;
        return (space);
    }
}