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
#' @param random_seed Seed for replicability.
#' @param max_exclusions maximum number of cases to be excluded
#'
#' @return Vector of row indeces to be excluded
#'
#' @examples
#' coefficient_computation <- function(data){
#' statistic <- cor(data$Sepal.Length, data$Petal.Width)
#' return(statistic)}
#'
#' filter <- loop_break(data = iris, statistic_computation = coefficient_computation, goal_value = 0.2, max_exclusions = 2)
#' print(filter)
#'
#' @export



# loop function -----------------------------------------------------------

loop_break = function(data = NA, goal_value = NA, statistic_computation = NA, max_exclusions = 3, random_seed = 42){
  set.seed(random_seed)
  indeces = seq(nrow(data)) #get row indeces (for making combinations)
  best_so_far = 42 #set random initial result
  mem_issue = FALSE
  first = TRUE #indicate it is the first round
  baseline_value = statistic_computation(data) #what is the full sample statistic
  assert("Goal value is the same as the full sample result", {goal_value != baseline_value}) #check that goal value is different
  minimize = ifelse(goal_value < baseline_value, TRUE, FALSE) #check whether to minimize or maximize
  for(i in 1:max_exclusions){ #first check all one-exclusion-possibilities, then all two-exclusion possibilities...
    indeces_combis = combn(indeces,i) #all row combinations

    for(combi in 1:ncol(indeces_combis)){ #for each combination of i rows
      current_exclusions = indeces_combis[,combi] #note the excluded rows
      current_result = statistic_computation(data[-current_exclusions,]) #compute the result for the sample without the rows

      if(first){ best_so_far = current_result; first = FALSE; best_exclusions = current_exclusions} #if it is the first computation, set result as best result

      if((minimize & current_result < best_so_far) | (!minimize & current_result > best_so_far) ) {#if result is better, update best result
        best_so_far = current_result
        best_exclusions = current_exclusions
      }
      if(combi %% 1000 == 0){cat('Step: ',i, '/', max_exclusions, ', Progress current step: ', combi, '/', ncol(indeces_combis), ', Target statistic: ', best_so_far, '\n', sep = '')}

    }
    if((minimize & best_so_far < goal_value) | (!minimize & best_so_far >  goal_value) ){#if goal value was reached; return best sample
      print('solution found'); return(sort(best_exclusions))
    }
    if(i == max_exclusions){print('no solution found'); return(NA)} #if no exclusions reach the goal value, return na
  }

}
