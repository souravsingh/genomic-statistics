data("Maupin")
names(maupin)

geneSet<- maupin$sig$EntrezID    
up_sig<- maupin$sig[maupin$sig$upDown == "up",]
d_sig<- maupin$sig[maupin$sig$upDown == "down",]
u_geneSet<- up_sig$EntrezID   
d_geneSet<- d_sig$EntrezID

enrichment_scores <- gsva(maupin$data, list(up = u_geneSet, down= d_geneSet), mx.diff=1,
               verbose=TRUE, abs.ranking=FALSE, is.gset.list.up.down=TRUE, parallel.sz = 1 )$es.obs
               
es.dif.ssg <- gsva(maupin, list(up = u_geneSet, down= d_geneSet),
                                                        verbose=TRUE, abs.ranking=FALSE, is.gset.list.up.down=TRUE,
                                                        method = "ssgsea")
                                                        
hist(enrichment_scores, main = "enrichment scores", xlab="es")
lines(density(enrichment_scores[,1:3]), col = "blue") # control samples
lines(density(enrichment_scores[,4:6]), col = "red") # TGFb samples
legend("topleft", c("Control","TGFb"), lty = 1, col=c("blue","red"), cex = 0.6)
