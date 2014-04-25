pp <- function (x = seq(0.5, 10, length.out=100), y = data.frame(sin(x), cos(x)), 
		  xlim=c(-1.1, 10.5), ylim=c(-1.1,10.5)) 
{
    plot(x, y[, 1], type = "l", col = "red", xlim = xlim, ylim = ylim)
    lines(x, y[, 2], col = "blue")
    lines(y[, 2], x, col = "green")
    lines(y[, 1], x, col = "yellow")
}
