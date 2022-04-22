tau_n_list=c(0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

rowsd <- function(x){
  return(sapply(1:dim(x)[1],function(i) sd(x[i,],na.rm = T)))
}

 #####################  plot triangle sbm ################################

load("EW-MB-triangle-tau-sbm-15taus.RData")
pdf("sbm-triangles.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(tau_n_list,rowMeans(res$sup_ew),type='l',ylab='CDF Error',xlab=expression(rho),cex.axis=0.75)
points(tau_n_list,rowMeans(res$sup_ew),pch=20,cex=0.3)
points(tau_n_list,rowMeans(res$sup_ew)-rowsd(res$sup_ew),pch='-',cex=0.5,cex.axis=0.7)
points(tau_n_list,rowMeans(res$sup_ew)+rowsd(res$sup_ew),pch='-',cex=0.5,cex.axis=0.7)

lines(tau_n_list,rowMeans(res$sup_add),type='l',col='blue')
points(tau_n_list,rowMeans(res$sup_add),pch=20,cex=0.3,col='blue')
points(tau_n_list,rowMeans(res$sup_add)-rowsd(res$sup_add),pch='-',cex=0.5,cex.axis=0.7,col='blue')
points(tau_n_list,rowMeans(res$sup_add)+rowsd(res$sup_add),pch='-',cex=0.5,cex.axis=0.7,col='blue')

lines(tau_n_list,rowMeans(res$sup_linqua),type='l',col='red')
points(tau_n_list,rowMeans(res$sup_linqua),pch=20,cex=0.3,col='red')
points(tau_n_list,rowMeans(res$sup_linqua)-rowsd(res$sup_linqua),pch='-',cex=0.5,cex.axis=0.7,col='red')
points(tau_n_list,rowMeans(res$sup_linqua)+rowsd(res$sup_linqua),pch='-',cex=0.5,cex.axis=0.7,col='red')

lines(tau_n_list,rowMeans(res$sup_gauprod),type='l',col='green')
points(tau_n_list,rowMeans(res$sup_gauprod),pch=20,cex=0.3,col='green')
points(tau_n_list,rowMeans(res$sup_gauprod)-rowsd(res$sup_gauprod),pch='-',cex=0.5,cex.axis=0.7,col='green')
points(tau_n_list,rowMeans(res$sup_gauprod)+rowsd(res$sup_gauprod),pch='-',cex=0.5,cex.axis=0.7,col='green')


lines(tau_n_list,rowMeans(res$sup_sub),type='l',col=8)
points(tau_n_list,rowMeans(res$sup_sub),pch=20,cex=0.3,col=8)
points(tau_n_list,rowMeans(res$sup_sub)-rowsd(res$sup_sub),pch='-',cex=0.5,cex.axis=0.7,col=8)
points(tau_n_list,rowMeans(res$sup_sub)+rowsd(res$sup_sub),pch='-',cex=0.5,cex.axis=0.7,col=8)

lines(tau_n_list,rowMeans(res$sup_lz),type='l',col=5)
points(tau_n_list,rowMeans(res$sup_lz),pch=20,cex=0.3,col=5)
points(tau_n_list,rowMeans(res$sup_lz)-rowsd(res$sup_lz),pch='-',cex=0.5,cex.axis=0.7,col=5)
points(tau_n_list,rowMeans(res$sup_lz)+rowsd(res$sup_lz),pch='-',cex=0.5,cex.axis=0.7,col=5)

lines(tau_n_list,rowMeans(res$sup_eg),type='l',col=7)
points(tau_n_list,rowMeans(res$sup_eg),pch=20,cex=0.3,col=7)
points(tau_n_list,rowMeans(res$sup_eg)-rowsd(res$sup_eg),pch='-',cex=0.5,cex.axis=0.7,col=7)
points(tau_n_list,rowMeans(res$sup_eg)+rowsd(res$sup_eg),pch='-',cex=0.5,cex.axis=0.7,col=7)

lines(tau_n_list,rowMeans(res$sup_ut),type='l',col='orange')
points(tau_n_list,rowMeans(res$sup_ut),pch=20,cex=0.3,col='orange')
points(tau_n_list,rowMeans(res$sup_ut)-rowsd(res$sup_ut),pch='-',cex=0.5,cex.axis=0.7,col='orange')
points(tau_n_list,rowMeans(res$sup_ut)+rowsd(res$sup_ut),pch='-',cex=0.5,cex.axis=0.7,col='orange')

legend(0.5,0.2,legend = c("EW", "MB-L","MB-Q","MB-M"),col = c('black',"blue","red","green"), lty = 1, cex = 0.65, bty='n')
legend(0.7,0.2,legend = c("SS","LS","EG","MB-L-apx"),col = c(8,5,7,"orange"), lty = 1, cex = 0.65, bty='n')

dev.off()

#####################  plot twostar sbm ################################

load("EW-MB-twostar-tau-sbm-15taus.RData")
pdf("sbm-vstar-centerbytnhat.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(tau_n_list,rowMeans(res$sup_ew,na.rm = T),type='l',ylim=c(0,0.25),ylab='CDF Error',xlab=expression(rho),cex.axis=0.75)
points(tau_n_list,rowMeans(res$sup_ew,na.rm = T),pch=20,cex=0.3)
points(tau_n_list,rowMeans(res$sup_ew,na.rm = T)-rowsd(res$sup_ew),pch='-',cex=0.5,cex.axis=0.7)
points(tau_n_list,rowMeans(res$sup_ew,na.rm = T)+rowsd(res$sup_ew),pch='-',cex=0.5,cex.axis=0.7)

lines(tau_n_list,rowMeans(res$sup_add,na.rm = T),type='l',col='blue')
points(tau_n_list,rowMeans(res$sup_add,na.rm = T),pch=20,cex=0.3,col='blue')
points(tau_n_list,rowMeans(res$sup_add,na.rm = T)-rowsd(res$sup_add),pch='-',cex=0.5,cex.axis=0.7,col='blue')
points(tau_n_list,rowMeans(res$sup_add,na.rm = T)+rowsd(res$sup_add),pch='-',cex=0.5,cex.axis=0.7,col='blue')

lines(tau_n_list,rowMeans(res$sup_linqua,na.rm = T),type='l',col='red')
points(tau_n_list,rowMeans(res$sup_linqua,na.rm = T),pch=20,cex=0.3,col='red')
points(tau_n_list,rowMeans(res$sup_linqua,na.rm = T)-rowsd(res$sup_linqua),pch='-',cex=0.5,cex.axis=0.7,col='red')
points(tau_n_list,rowMeans(res$sup_linqua,na.rm = T)+rowsd(res$sup_linqua),pch='-',cex=0.5,cex.axis=0.7,col='red')

lines(tau_n_list,rowMeans(res$sup_gauprod,na.rm = T),type='l',col='green')
points(tau_n_list,rowMeans(res$sup_gauprod,na.rm = T),pch=20,cex=0.3,col='green')
points(tau_n_list,rowMeans(res$sup_gauprod,na.rm = T)-rowsd(res$sup_gauprod),pch='-',cex=0.5,cex.axis=0.7,col='green')
points(tau_n_list,rowMeans(res$sup_gauprod,na.rm = T)+rowsd(res$sup_gauprod),pch='-',cex=0.5,cex.axis=0.7,col='green')


lines(tau_n_list,rowMeans(res$sup_sub),type='l',col=8)
points(tau_n_list,rowMeans(res$sup_sub),pch=20,cex=0.3,col=8)
points(tau_n_list,rowMeans(res$sup_sub)-rowsd(res$sup_sub),pch='-',cex=0.5,cex.axis=0.7,col=8)
points(tau_n_list,rowMeans(res$sup_sub)+rowsd(res$sup_sub),pch='-',cex=0.5,cex.axis=0.7,col=8)

lines(tau_n_list,rowMeans(res$sup_lz),type='l',col=5)
points(tau_n_list,rowMeans(res$sup_lz),pch=20,cex=0.3,col=5)
points(tau_n_list,rowMeans(res$sup_lz)-rowsd(res$sup_lz),pch='-',cex=0.5,cex.axis=0.7,col=5)
points(tau_n_list,rowMeans(res$sup_lz)+rowsd(res$sup_lz),pch='-',cex=0.5,cex.axis=0.7,col=5)

lines(tau_n_list,rowMeans(res$sup_eg),type='l',col=7)
points(tau_n_list,rowMeans(res$sup_eg),pch=20,cex=0.3,col=7)
points(tau_n_list,rowMeans(res$sup_eg)-rowsd(res$sup_eg),pch='-',cex=0.5,cex.axis=0.7,col=7)
points(tau_n_list,rowMeans(res$sup_eg)+rowsd(res$sup_eg),pch='-',cex=0.5,cex.axis=0.7,col=7)

lines(tau_n_list,rowMeans(res$sup_ut),type='l',col='orange')
points(tau_n_list,rowMeans(res$sup_ut),pch=20,cex=0.3,col='orange')
points(tau_n_list,rowMeans(res$sup_ut)-rowsd(res$sup_ut),pch='-',cex=0.5,cex.axis=0.7,col='orange')
points(tau_n_list,rowMeans(res$sup_ut)+rowsd(res$sup_ut),pch='-',cex=0.5,cex.axis=0.7,col='orange')

legend(0.5,0.25,legend = c("EW", "MB-L","MB-Q","MB-M"),col = c('black',"blue","red","green"), lty = 1, cex = 0.75, bty='n')
legend(0.7,0.25,legend = c("SS","LS","EG","MB-L-apx"),col = c(8,5,7,"orange"), lty = 1, cex = 0.75, bty='n')

dev.off()


#####################  plot triangle  smg ################################
load("EW-MB-triangle-tau-smg-15taus.RData")
pdf("smg-triangles.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(tau_n_list,rowMeans(res$sup_ew,na.rm = T),type='l',ylab="CDF Error",xlab=expression(rho),cex.axis=0.75)
points(tau_n_list,rowMeans(res$sup_ew,na.rm = T),pch=20,cex=0.3)
points(tau_n_list,rowMeans(res$sup_ew,na.rm = T)-rowsd(res$sup_ew),pch='-',cex=0.5,cex.axis=0.7)
points(tau_n_list,rowMeans(res$sup_ew,na.rm = T)+rowsd(res$sup_ew),pch='-',cex=0.5,cex.axis=0.7)

lines(tau_n_list,rowMeans(res$sup_add,na.rm = T),type='l',col='blue')
points(tau_n_list,rowMeans(res$sup_add,na.rm = T),pch=20,cex=0.3,col='blue')
points(tau_n_list,rowMeans(res$sup_add,na.rm = T)-rowsd(res$sup_add),pch='-',cex=0.5,cex.axis=0.7,col='blue')
points(tau_n_list,rowMeans(res$sup_add,na.rm = T)+rowsd(res$sup_add),pch='-',cex=0.5,cex.axis=0.7,col='blue')

lines(tau_n_list,rowMeans(res$sup_linqua,na.rm = T),type='l',col='red')
points(tau_n_list,rowMeans(res$sup_linqua,na.rm = T),pch=20,cex=0.3,col='red')
points(tau_n_list,rowMeans(res$sup_linqua,na.rm = T)-rowsd(res$sup_linqua),pch='-',cex=0.5,cex.axis=0.7,col='red')
points(tau_n_list,rowMeans(res$sup_linqua,na.rm = T)+rowsd(res$sup_linqua),pch='-',cex=0.5,cex.axis=0.7,col='red')

lines(tau_n_list,rowMeans(res$sup_gauprod,na.rm = T),type='l',col='green')
points(tau_n_list,rowMeans(res$sup_gauprod,na.rm = T),pch=20,cex=0.3,col='green')
points(tau_n_list,rowMeans(res$sup_gauprod,na.rm = T)-rowsd(res$sup_gauprod),pch='-',cex=0.5,cex.axis=0.7,col='green')
points(tau_n_list,rowMeans(res$sup_gauprod,na.rm = T)+rowsd(res$sup_gauprod),pch='-',cex=0.5,cex.axis=0.7,col='green')


lines(tau_n_list,rowMeans(res$sup_sub,na.rm = T),type='l',col=8)
points(tau_n_list,rowMeans(res$sup_sub,na.rm = T),pch=20,cex=0.3,col=8)
points(tau_n_list,rowMeans(res$sup_sub,na.rm = T)-rowsd(res$sup_sub),pch='-',cex=0.5,cex.axis=0.7,col=8)
points(tau_n_list,rowMeans(res$sup_sub,na.rm = T)+rowsd(res$sup_sub),pch='-',cex=0.5,cex.axis=0.7,col=8)

lines(tau_n_list,rowMeans(res$sup_lz,na.rm = T),type='l',col=5)
points(tau_n_list,rowMeans(res$sup_lz,na.rm = T),pch=20,cex=0.3,col=5)
points(tau_n_list,rowMeans(res$sup_lz,na.rm = T)-rowsd(res$sup_lz),pch='-',cex=0.5,cex.axis=0.7,col=5)
points(tau_n_list,rowMeans(res$sup_lz,na.rm = T)+rowsd(res$sup_lz),pch='-',cex=0.5,cex.axis=0.7,col=5)

lines(tau_n_list,rowMeans(res$sup_eg,na.rm = T),type='l',col=7)
points(tau_n_list,rowMeans(res$sup_eg,na.rm = T),pch=20,cex=0.3,col=7)
points(tau_n_list,rowMeans(res$sup_eg,na.rm = T)-rowsd(res$sup_eg),pch='-',cex=0.5,cex.axis=0.7,col=7)
points(tau_n_list,rowMeans(res$sup_eg,na.rm = T)+rowsd(res$sup_eg),pch='-',cex=0.5,cex.axis=0.7,col=7)

lines(tau_n_list,rowMeans(res$sup_ut,na.rm = T),type='l',col='orange')
points(tau_n_list,rowMeans(res$sup_ut,na.rm = T),pch=20,cex=0.3,col='orange')
points(tau_n_list,rowMeans(res$sup_ut,na.rm = T)-rowsd(res$sup_ut),pch='-',cex=0.5,cex.axis=0.7,col='orange')
points(tau_n_list,rowMeans(res$sup_ut,na.rm = T)+rowsd(res$sup_ut),pch='-',cex=0.5,cex.axis=0.7,col='orange')

legend(0.5,0.25,legend = c("EW", "MB-L","MB-Q","MB-M"),col = c('black',"blue","red","green"), lty = 1, cex = 0.65, bty='n')
legend(0.7,0.25,legend = c("SS","LS","EG","MB-L-apx"),col = c(8,5,7,"orange"), lty = 1, cex = 0.65, bty='n')

dev.off()

#####################  plot twostar smg ################################

load("EW-MB-twostar-tau-smg-15taus.RData")
pdf("smg-vstar.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(tau_n_list,rowMeans(res$sup_ew),type='l',ylim=c(0,0.3),ylab='CDF Error',xlab=expression(rho),cex.axis=0.75)
points(tau_n_list,rowMeans(res$sup_ew),pch=20,cex=0.3)
points(tau_n_list,rowMeans(res$sup_ew)-rowsd(res$sup_ew),pch='-',cex=0.5,cex.axis=0.7)
points(tau_n_list,rowMeans(res$sup_ew)+rowsd(res$sup_ew),pch='-',cex=0.5,cex.axis=0.7)

lines(tau_n_list,rowMeans(res$sup_add),type='l',col='blue')
points(tau_n_list,rowMeans(res$sup_add),pch=20,cex=0.3,col='blue')
points(tau_n_list,rowMeans(res$sup_add)-rowsd(res$sup_add),pch='-',cex=0.5,cex.axis=0.7,col='blue')
points(tau_n_list,rowMeans(res$sup_add)+rowsd(res$sup_add),pch='-',cex=0.5,cex.axis=0.7,col='blue')

lines(tau_n_list,rowMeans(res$sup_linqua),type='l',col='red')
points(tau_n_list,rowMeans(res$sup_linqua),pch=20,cex=0.3,col='red')
points(tau_n_list,rowMeans(res$sup_linqua)-rowsd(res$sup_linqua),pch='-',cex=0.5,cex.axis=0.7,col='red')
points(tau_n_list,rowMeans(res$sup_linqua)+rowsd(res$sup_linqua),pch='-',cex=0.5,cex.axis=0.7,col='red')

lines(tau_n_list,rowMeans(res$sup_gauprod),type='l',col='green')
points(tau_n_list,rowMeans(res$sup_gauprod),pch=20,cex=0.3,col='green')
points(tau_n_list,rowMeans(res$sup_gauprod)-rowsd(res$sup_gauprod),pch='-',cex=0.5,cex.axis=0.7,col='green')
points(tau_n_list,rowMeans(res$sup_gauprod)+rowsd(res$sup_gauprod),pch='-',cex=0.5,cex.axis=0.7,col='green')


lines(tau_n_list,rowMeans(res$sup_sub),type='l',col=8)
points(tau_n_list,rowMeans(res$sup_sub),pch=20,cex=0.3,col=8)
points(tau_n_list,rowMeans(res$sup_sub)-rowsd(res$sup_sub),pch='-',cex=0.5,cex.axis=0.7,col=8)
points(tau_n_list,rowMeans(res$sup_sub)+rowsd(res$sup_sub),pch='-',cex=0.5,cex.axis=0.7,col=8)

lines(tau_n_list,rowMeans(res$sup_lz),type='l',col=5)
points(tau_n_list,rowMeans(res$sup_lz),pch=20,cex=0.3,col=5)
points(tau_n_list,rowMeans(res$sup_lz)-rowsd(res$sup_lz),pch='-',cex=0.5,cex.axis=0.7,col=5)
points(tau_n_list,rowMeans(res$sup_lz)+rowsd(res$sup_lz),pch='-',cex=0.5,cex.axis=0.7,col=5)

lines(tau_n_list,rowMeans(res$sup_eg),type='l',col=7)
points(tau_n_list,rowMeans(res$sup_eg),pch=20,cex=0.3,col=7)
points(tau_n_list,rowMeans(res$sup_eg)-rowsd(res$sup_eg),pch='-',cex=0.5,cex.axis=0.7,col=7)
points(tau_n_list,rowMeans(res$sup_eg)+rowsd(res$sup_eg),pch='-',cex=0.5,cex.axis=0.7,col=7)

lines(tau_n_list,rowMeans(res$sup_ut),type='l',col='orange')
points(tau_n_list,rowMeans(res$sup_ut),pch=20,cex=0.3,col='orange')
points(tau_n_list,rowMeans(res$sup_ut)-rowsd(res$sup_ut),pch='-',cex=0.5,cex.axis=0.7,col='orange')
points(tau_n_list,rowMeans(res$sup_ut)+rowsd(res$sup_ut),pch='-',cex=0.5,cex.axis=0.7,col='orange')

legend(0.5,0.3,legend = c("EW", "MB-L","MB-Q","MB-M"),col = c('black',"blue","red","green"), lty = 1, cex = 0.65, bty='n')
legend(0.7,0.3,legend = c("SS","LS","EG","MB-L-apx"),col = c(8,5,7,"orange"), lty = 1, cex = 0.65, bty='n')


dev.off()

