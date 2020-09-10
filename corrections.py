# this code goes into code cell 1.1

def random_optimize(domain,  fitness_function, n=9999):    
    best_score = None
    best_sol = None
    for iter in range(n):
        # Create a random solution
        # domain[i] represents a range of all valid choices for position i
        r=[random.choice(domain[i])
                      for i in range(len(domain))]
    
        # Get the cost
        cost=fitness_function(r)
    
        # Compare it to the best one so far
        if best_score is None or cost < best_score:
            best_score = cost
            best_sol = r

    return (best_sol, best_score, n)
    
    
# this code goes into code cell 1.2
def hillclimb_optimize(domain, fitness_function):
    # Create an initial random solution
    best_sol=[random.choice(domain[i])
              for i in range(len(domain))]

    # Main loop: until no better neighbors
    iterations = 0
    while 1:
        best_score=fitness_function(best_sol)        
        current_score = best_score
        
        iterations += 1
        
        # Create list of all neighboring solutions
        neighbors=[]

        for j in range(len(domain)):
            # One away in each direction
            # if current best-SOl[j] is not the smallest in the domain 
            # - we can subtract 1
            if best_sol[j]>domain[j][0]: 
                neighbors.append(best_sol[0:j]+[best_sol[j]-1]+best_sol[j+1:])
            # if current best-SOl[j] is not the largest in the domain 
            # - we can add 1 
            if best_sol[j]< domain[j][1] - 1:
                neighbors.append(best_sol[0:j]+[best_sol[j]+1]+best_sol[j+1:])

        # See what the best solution amongst the neighbors is        
        for j in range(len(neighbors)):
            cost = fitness_function(neighbors[j])
            if cost < best_score:
                best_score = cost
                best_sol = neighbors[j]

        # If there's no improvement, then we've reached the local min
        if best_score == current_score:
            break
            
        if best_score == 0:
            break

    return (best_sol, best_score, iterations)
    

# this code goes into code cell 1.3
def annealing_optimize(domain, fitness_function,
                      T=10000.0,cool=0.95,step=1):
    # Create an initial random solution
    best_sol =[random.choice(domain[i])
              for i in range(len(domain))]
    
    iterations = 0
    while T > 0.1:
        iterations += 1
        
        # Choose one of the indices at random
        i = random.randint(0,len(domain)-1) #randint selects INCLUDING the max

        # Choose a random direction to change it
        dir = random.randint(-step,step)

        # Create a new solution with one of the values changed into a random direction
        sol_b = best_sol[:]
        sol_b [i] += dir
        
        # check that it is within allowed domain
        # if sol_b[i] smaller than the min allowed value
        # replace it with the smallest from the domain
        if sol_b [i] < domain[i][0]:
            sol_b [i] = domain[i][0]
        # if sol_b[i] larger than the max allowed value
        # replace it with the largest from the domain
        elif sol_b [i] > domain[i][1]-1:
            sol_b [i] = domain[i][1] - 1

        # Calculate the cost of the current_best solution
        best_score = fitness_function(best_sol)
        
        if best_score == 0:
            break
            
        # Calculate score of a new solution   
        b_score = fitness_function(sol_b)
        
        # If T is really large, exp -> 0, and p -> 1
        # If T is small, exp is defined by difference between new and old score
        # If T is very small, exp -> infinity, p -> 0 
        exp = -(b_score - best_score)/T        
        p = pow (math.e, exp)

        # Is it better, 
        # or does it make the probability cutoff?
        if b_score < best_score or p > random.random():
            best_sol = sol_b        
        
        # Decrease the temperature
        T = T*cool

    return (best_sol, fitness_function(best_sol), iterations)
 
 
# this code goes into code cell 1.4  
def genetic_optimize (domain, fitness_function,
                    popsize=100, step=1,
                    mut_prob=0.2, elite_ratio=0.2, n=100):
    # Mutation Operation
    def mutate(vec):
        # choose random index in the solution
        i=random.randint(0,len(domain)-1)
        # Rooll the dice and with equal probability
        # either increase or decrease value at i by step
        if random.random()<0.5 and (vec[i] - step) >= domain[i][0]:
            return vec[0:i]+[vec[i]-step]+vec[i+1:]
        elif (vec[i]+step) < domain[i][1]:
            return vec[0:i]+[vec[i]+step]+vec[i+1:]
        return vec

    # Crossover Operation
    def crossover(r1,r2):
        i=random.randint(0,len(domain)-1)
        return r1[0:i]+r2[i:]

    # Build the initial population of random solutions
    pop=[]
    for i in range(popsize):
        vec = [random.choice(domain[i])
               for i in range(len(domain))]
        pop.append(vec)

    # How many winners from each generation?
    topelite=int(elite_ratio*popsize)

    # Main loop
    generations = 0
    for i in range(n):
        # This is the list of all solutions in the population together with their fitness scores
        scores = [(fitness_function(sol_vect), sol_vect) for sol_vect in pop]
        
        # We sort this list by score
        scores.sort()

        # See what is current top best score
        if scores[0][0] == 0:
            break
        
        generations += 1
        
        # this is a list of just the solutions extracted from the sorted by score
        ranked_solutions = [sol_vect for (cost, sol_vect) in scores]

        # Build next gen population
        # Start with the pure winners
        pop = ranked_solutions[0:topelite]

        # Add mutated and bred forms of the winners
        while len(pop) < popsize:
            if random.random() < mut_prob:
                # Mutation. Select random individual from elite group and mutate
                c = random.randint(0,topelite-1)
                pop.append(mutate(ranked_solutions[c]))
            else:
                # Crossover. Select 2 random individuals from elite group and cross
                c1 = random.randint(0, topelite-1)
                c2 = random.randint(0, topelite-1)
                pop.append(crossover(ranked_solutions[c1],ranked_solutions[c2]))
                
    # After n generations return the top scored solution from the sorted list of scores
    return (scores[0][1], scores[0][0], generations) 