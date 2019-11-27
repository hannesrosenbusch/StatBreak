#' Find smallest subset to exclude from sample for effect/pattern to disappear
#'
#' The function iteratively learns which observations should at least be excluded from
#' the data to reach a conservative 'goal value' for the statistic of interest.
#' It does so by relying on a genetic algorithm, which efficiently explores the (usually vast)
#' space of possible subsets. The result can uncover impactful subsamples and fuel discussions of robustness.
#' Necessary arguments include the dataframe,
#' a function to compute the statistic of interest ('statistic_computation' see examples),
#' and the goal value of interest.
#'
#' @param data A data.frame containing the observations as rows.
#' @param goal_value This conservative value (e.g., small effect size) is targeted.
#' @param statistic_computation A formula which has 'data' as input and returns the statistic of interest.
#' @param pop Number of 'individuals' in each generation of the genetic algorithm.
#' @param max_generations Maximum number of generations that the algorithm generates.
#' @param exclusion_cost Used to calibrate fitness function.
#' @param prop_included_cases Initial proportion of included cases (e.g. .90).
#' @param chance_of_mutation Chance that a gene mutates, higher is slower but more accurate (e.g. .02).
#' @param stop_search After how many generations without change is the 'converged' result returned.
#' @param random_seed Seed for replicability.
#' @param max_exclusions maximum number of cases to be excluded
#'
#' @return Vector of row indeces of rows to be excluded
#'
#' @examples
#' coefficient_computation <- function(data){
#' statistic <- cor(data$Sepal.Length, data$Petal.Width)
#' return(statistic)}
#'
#' filter <- stat_break(data = iris, statistic_computation = coefficient_computation, goal_value = 0.2, max_generations = 500)
#' print(filter)
#'
#' @export



# break function -----------------------------------------------------------

stat_break = function(data = NULL,
                   goal_value = NULL, #how many cases need to be excluded to reach this value
                   statistic_computation = NULL, #formula for the computation of the statistic
                   max_exclusions = NULL,

                   pop = 1000, #population size of each generation (number datasets), higher is slower but more accurate min 500
                   max_generations = 2000, #maximum number of generations. higher is slower but more accurate min 100
                   exclusion_cost = 0.01, #stable exclusion costs, should not need tuning
                   prop_included_cases = 0.95, #initial proportion of included cases (0-1), lower is much slower but more accurate
                   chance_of_mutation = 0.02, #chance that a gene mutates, higher is slower but more accurate max 0.1
                   stop_search = 200,#after how many generations without improvements is result returned
                   random_seed = 42){


  # setup -------------------------------------------------------------------
  set.seed(random_seed)
  if(is.null(max_exclusions)){
  max_exclusions = floor(nrow(data)*0.5)
  print(paste("max_exclusions argument was not provided. Using ", max_exclusions, " as limit."))
  }
  if(is.null(data)){
    stop("Please provide the data under the 'data' argument.")
  }
  if(is.null(statistic_computation)){
    stop("A function must be provided as the 'statistic_computation' argument which takes data as input and returns the statistic of interest.")
  }
  if(is.null(goal_value)){
    stop("A numeric value must be provided as the 'goal_value' argument which indcates the conservative point of interest for the researcher.")
  }

  # define custom version of genalg function ---------------------------------------
  #original version here: https://github.com/cran/genalg/tree/master/R
  rbga.bin <- function(
    size=NA,
    suggestions=NULL,
    popSize=NA, iters=NA,
    mutationChance=NA,
    elitism=NA, zeroToOneRatio=NA,
    monitorFunc=NULL, evalFunc=NULL,
    showSettings=FALSE, verbose=FALSE,
    parentProb= NA, stop_search = NA
  ) {
    if (is.null(evalFunc)) {
      stop("An evaluation function must be provided. See the evalFunc parameter.");
    }

    vars = size;
    if (is.na(mutationChance)) {
      mutationChance = 1/(vars+1);
    }
    if (is.na(elitism)) {
      elitism = floor(popSize/5)
    }

    if (verbose) cat("Testing the sanity of parameters...\n");
    if (popSize < 5) {
      stop("The population size must be at least 5.");
    }
    if (iters < 1) {
      stop("The number of iterations must be at least 1.");
    }
    if (!(elitism < popSize)) {
      stop("The population size must be greater than the elitism.");
    }

    if (showSettings) {
      if (verbose) cat("The start conditions:\n");
      result = list(size=size, suggestions=suggestions,
                    popSize=popSize, iters=iters,
                    elitism=elitism, mutationChance=mutationChance);
      class(result) = "rbga";

      cat(summary(result));
    } else {
      if (verbose) cat("Not showing GA settings...\n");
    }

    if (vars > 0) {
      if (!is.null(suggestions)) {
        if (verbose) cat("Adding suggestions to first population...\n");
        population = matrix(nrow=popSize, ncol=vars);
        suggestionCount = dim(suggestions)[1]
        for (i in 1:suggestionCount) {
          population[i,] = suggestions[i,]
        }
        if (verbose) cat("Filling others with random values in the given domains...\n");
        for (child in (suggestionCount+1):popSize) {
          population[child,] = sample(c(rep(0,zeroToOneRatio),1), vars, replace=TRUE);
          while (sum(population[child,]) == 0) {
            population[child,] = sample(c(rep(0,zeroToOneRatio),1), vars, replace=TRUE);
          }
        }
      } else {
        if (verbose) cat("Starting with random values in the given domains...\n");
        # start with an random population
        population = matrix(nrow=popSize, ncol=vars);
        # fill values
        for (child in 1:popSize) {
          population[child,] = sample(c(rep(0,zeroToOneRatio),1), vars, replace=TRUE);
          while (sum(population[child,]) == 0) {
            population[child,] = sample(c(rep(0,zeroToOneRatio),1), vars, replace=TRUE);
          }
        }
      }

      # do iterations
      stability = 1
      bestsolution = 999999
      bestEvals = rep(NA, iters);
      meanEvals = rep(NA, iters);
      evalVals = rep(NA, popSize);
      for (iter in 1:iters) {

        if (verbose) cat(paste("Starting iteration", iter, "\n"));

        # calculate each object
        if (verbose) cat("Calucating evaluation values... ");
        for (object in 1:popSize) {
          if (is.na(evalVals[object])) {
            evalVals[object] = evalFunc(population[object,]);
            if (verbose) cat(".");
          }
        }


        if(stability == stop_search){
          #print(paste(stop_search,"iterations without change"))
          result = list(type="binary chromosome", size=size,
                        popSize=popSize, iters=iters, suggestions=suggestions,
                        population=population, elitism=elitism, mutationChance=mutationChance,
                        evaluations=evalVals, best=bestEvals, mean=meanEvals, stability = stability, stop_search = stop_search);
          class(result) = "rbga";

          return(result);
        }


        bestEvals[iter] = min(evalVals);

        if (min(evalVals) == bestsolution){
          stability = stability +1;
          bestsolution = min(evalVals);
        }else{
          bestsolution = min(evalVals);
          stability = 1
        }


        meanEvals[iter] = mean(evalVals);
        if (verbose) cat(" done.\n");

        if (!is.null(monitorFunc)) {
          if (verbose) cat("Sending current state to rgba.monitor()...\n");
          # report on GA settings
          result = list(type="binary chromosome", size=size,
                        popSize=popSize, iter=iter, iters=iters,
                        population=population, elitism=elitism, mutationChance=mutationChance,
                        evaluations=evalVals, best=bestEvals, mean=meanEvals);
          class(result) = "rbga";

          monitorFunc(result);
        }

        if (iter < iters) { # ok, must create the next generation
          if (verbose) cat("Creating next generation...\n");
          newPopulation = matrix(nrow=popSize, ncol=vars);
          newEvalVals = rep(NA, popSize);

          if (verbose) cat("  sorting results...\n");
          sortedEvaluations = sort(evalVals, index=TRUE);
          sortedPopulation  = matrix(population[sortedEvaluations$ix,], ncol=vars);

          # save the best
          if (elitism > 0) {
            if (verbose) cat("  applying elitism...\n");
            newPopulation[1:elitism,] = sortedPopulation[1:elitism,];
            newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
          } # ok, save nothing

          # fill the rest by doing crossover
          if (vars > 1) {
            if (verbose) cat("  applying crossover...\n");
            for (child in (elitism+1):popSize) {
              # ok, pick two random parents
              parentIDs = sample(1:popSize, 2, prob= parentProb)
              parents = sortedPopulation[parentIDs,];
              crossOverPoint = sample(0:vars,1);
              if (crossOverPoint == 0) {
                newPopulation[child, ] = parents[2,]
                newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
              } else if (crossOverPoint == vars) {
                newPopulation[child, ] = parents[1,]
                newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
              } else {
                newPopulation[child, ] =
                  c(parents[1,][1:crossOverPoint],
                    parents[2,][(crossOverPoint+1):vars])
                while (sum(newPopulation[child,]) == 0) {
                  newPopulation[child,] = sample(c(rep(0,zeroToOneRatio),1), vars, replace=TRUE);
                }
              }
            }
          } else { # otherwise nothing to crossover
            if (verbose) cat("  cannot crossover (#vars=1), using new randoms...\n");
            # random fill the rest
            newPopulation[(elitism+1):popSize,] =
              sortedPopulation[sample(1:popSize, popSize-elitism),];
          }

          population = newPopulation;
          evalVals   = newEvalVals;

          # do mutation
          if (mutationChance > 0) {
            if (verbose) cat("  applying mutations... ");
            mutationCount = 0;
            for (object in (elitism+1):popSize) {
              for (var in 1:vars) {
                if (runif(1) < mutationChance) { # ok, do mutation

                  population[object,var] = sample(c(rep(0,zeroToOneRatio),1), 1);
                  mutationCount = mutationCount + 1;
                }
              }
            }

            if (verbose) cat(paste(mutationCount, "mutations applied\n"));
          }
        }


      }
    }

    # report on GA settings
    result = list(type="binary chromosome", size=size,
                  popSize=popSize, iters=iters, suggestions=suggestions,
                  population=population, elitism=elitism, mutationChance=mutationChance,
                  evaluations=evalVals, best=bestEvals, mean=meanEvals);
    class(result) = "rbga";

    return(result);
  }


  # check if maximize or minimize -------------------------------------------
  v = statistic_computation(data)

  if(v > goal_value){
    maximize = FALSE
  }else if(v < goal_value){#not optimal to always recompute, should also check if numeric
    maximize = TRUE
  }else if(v == goal_value){
    stop("Your goal value is equal to the observed statistic in the full sample")
  }


  # implement evaluate function based on statistic_computation ---------------
  if(maximize){
  evaluate = function(string=c()){
        if(sum(string)> max_exclusions){
          return(999999)
        }
        x = statistic_computation(data[0==string,])
        return(sum(string)/length(string)*exclusion_cost + 1/min(x, exclusion_cost^4 * x + goal_value))
      }
  } else if(!maximize){

    evaluate = function(string=c()){
      if(sum(string)> max_exclusions){
        return(999999)
      }
      x = statistic_computation(data[0==string,])
      return(sum(string)/length(string)*exclusion_cost + max(x, exclusion_cost^4 * x + goal_value))
       }
  }


  # implement monitor function based on statistic_computation ---------------
  monitor = function(alg){
    minEval = min(alg$evaluations)
    filter = alg$evaluations == minEval
    bestObjectCount = sum(rep(1, alg$popSize)[filter])
    if (bestObjectCount > 1) {
      bestSolution = alg$population[filter,][1,]
    } else {
      bestSolution = alg$population[filter,]}
    x = statistic_computation(data[0==bestSolution,])
    cat('Generations w.o. change: ', alg$stability, '/', alg$stop_search, ',', 'dropped rows: ', 'sum(bestSolution)', ',', 'target statistic: ', x)
    # print(alg$iter)
    # print(paste("dropped observations:", sum(bestSolution)))
    # print(paste("observed statistic:", x))
  }



  # function for returning best solution after alg finished ---------------------------------------------
  best <- function(alg) {
    minEval = min(alg$evaluations, na.rm = T)
    filter = alg$evaluations == minEval
    bestObjectCount = sum(rep(1, alg$popSize)[filter])
    # ok, deal with the situation that more than one object is best
    if (bestObjectCount > 1) {
      bestSolution = alg$population[filter,][1,]
    } else {
      bestSolution = alg$population[filter,]}
    bestSolution}

  # run genetic alg with binary genes ------------------------------------------------------------
  alg = rbga.bin(size=nrow(data),
                 suggestions=NULL,
                 popSize=pop, iters=max_generations,
                 mutationChance=chance_of_mutation,
                 elitism=floor(0.1*pop), zeroToOneRatio=(prop_included_cases *nrow(data)) / ((1-prop_included_cases)* nrow(data) + 1^-10),#genes,
                 monitorFunc=monitor, evalFunc=evaluate,
                 showSettings=FALSE, verbose=FALSE, stop_search = stop_search,
                 parentProb= dnorm(1:pop, mean=0, sd=(pop/3)))

  # output best filter -------------------------------------------------------------
  solution = best(alg)
  print('Exclude the following observations (rows) for a less interesting finding:')
  print(which(solution == 1))
  return(sort(which(solution == 1))) #returns row indeces
}
