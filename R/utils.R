#####################################################
# Name: utils.R                                     #
# Author: Chris Jewell <c.jewell@lancaster.ac.uk>   #
# Created: 20161206                                 #
# Copyright: Chris Jewell 2016                      #
# Purpose: Miscellaneous helper functions           #
#####################################################


isFiniteInteger = function(a)
{
  return(is.finite(a) & isTRUE(all.equal(a, as.integer(a))))
}
isFiniteLogical = function(a)
{
  return(is.finite(a) & is.logical(a))
}
isFinitePositive = function(a)
{
  return(is.finite(a) & a > 0)
}
