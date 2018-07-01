package strategies;

import automail.StorageTube;

public class MyRobotBehaviour implements IRobotBehaviour {
    /***
     * Robot is a dummy - no abortion of mail runs
     */
    
    public MyRobotBehaviour(boolean strong) {
    }

    @Override
    public void startDelivery() {
    }

    @Override
    public boolean returnToMailRoom(StorageTube tube) {
        return false;
    }

    @Override
    public void priorityArrival(int priority, int weight) {
    }

}
