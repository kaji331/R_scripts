primest <- function(n){
    p <- 2:n
    i <- 1
    while (p[i] <= sqrt(n)) {
        p <-  p[p %% p[i] != 0 | p==p[i]]
        i <- i+1
    }
    p
}
