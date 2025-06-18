MatingModel <- R6::R6Class(
    "MatingModel",
    public = list(
        initialize = function(type, pairs = NULL) {
            checkmate::assert_choice(type, choices = c("random", "assortative"))
            checkmate::assert_class(pairs, "matrix")

            private$type <- type
            private$pairs <- pairs
        },

        get_type = function() return(private$type),
        get_pairs = function() return(private$pairs)
    ),
    private = list(
        type = NULL,
        pairs = NULL
    )
)