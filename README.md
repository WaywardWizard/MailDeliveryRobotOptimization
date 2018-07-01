# About this project
Project code responsible for the simulation is not commited here - the only files we were able to edit for this project were the two in this repository. This is here to demonstrate a basic competency with Java

## The project assumes the scenario;
- We have a high rise residential building, with a floor which recieves all residents mail.
- The building has two delivery robots that have a limited carry capacity (4 items) and limited maximum item weight that they may carry. 
- All mail has a weight and priority
- We are given some delivery costing function, f(priority, time taken to deliever)
- Mail arrives at random times with normally distributed weights about some mean

## The problem this code solves is;
- *Define robot pair behaviour that minimizes the given cost function*

## Roughly, whats the strategy?
The approach taken is to brute force carry set combinations for each robot. Known information about each robot and the world state is used to reduce the search space. Some (stated) assumptions are made to further reduce the search space. 

The project brief is included in pdf format for further detail
