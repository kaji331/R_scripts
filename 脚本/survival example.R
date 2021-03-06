# example for survival analysis
library(survival)
time <- c(1, 5, 10, 15, 30, 3, 5, 7, 9, 18, 12, 19, 20, 21, 33)
status <- c(0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0)
group <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
plot(survfit(Surv(time[1:5],status[1:5])~group[1:5]), col = "dark green", xlim = c(0, 33), conf.int = F, main = "K-M estimate for data", xlab = "Time", ylab = "Probability of survival")
lines(survfit(Surv(time[6:10],status[6:10])~group[6:10]), col = "dark red")
lines(survfit(Surv(time[11:15],status[11:15])~group[11:15]), col = "dark blue")
survdiff(Surv(time[c(1:5,11:15)],status[c(1:5,11:15)])~group[c(1:5,11:15)])
survdiff(Surv(time[6:15],status[6:15])~group[6:15], rho = 0)
survdiff(Surv(time[6:15],status[6:15])~group[6:15], rho = 1)
coxph(Surv(time[6:15],status[6:15])~group[6:15])
t <- function(x) {
	x + 1
}
print(t(9))

# ======
print(t(9))
print(t(9))
# ======
