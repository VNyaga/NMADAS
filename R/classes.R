#' @name nmadasdata-class
#' @title Class nmadasdata
#' @description A nmadasdata class in the NMADAS package.
#' @docType class
#' @slot data A data-frame with no missing values containg TP, TN, FP, FN, study and test names.
#' @slot S.ID Study identifier
#' @slot T.ID Test identifier
#' @slot labels A vector with test labels.
#' @family nmadas
#' @seealso \link{cdtamodel}
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}

setClass(Class="nmadasdata",
          representation=representation(
              data = 'data.frame',
              labels = 'character',
              S.ID = 'character',
              T.ID = 'character'))

#' @name nmadasmodel-class
#' @title Class nmadasmodel
#' @description A nmadasmodel class in the NMADAS package.
#' @docType class
#' @slot model Model specification in a character format.
#' @family nmadas
#' @seealso \link{nmamodel}
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}

setClass(Class="nmadasmodel",
         representation=representation(
           model = 'character'))

#' @name nmadasfit-class
#' @title Class nmadasfit
#' @description A nmadasfit class in the NMADAS package.
#' @docType class
#' @slot data data A data-frame with no missing values containg TP, TN, FP, FN, study and test names.
#' @slot comparator Name of comparator test.
#' @slot fit an object of class stanfit returned by the function sampling.
#' @slot S.ID Study identifier
#' @slot T.ID Test identifier
#' @family nmadas
#' @seealso \link{fit}
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
#' @importClassesFrom rstan stanmodel stanfit

setClass(Class="nmadasfit",
         representation=representation(
            data='data.frame',
            S.ID = 'character',
            T.ID = 'character',
            labels = 'character',
            comparator = 'character',
            fit = 'stanfit'))
