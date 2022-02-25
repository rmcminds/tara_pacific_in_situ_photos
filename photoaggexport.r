input_prefix <- '~/Dropbox/Tara/unsupervised_analyses'
output_prefix <- '~/Dropbox/Tara/unsupervised_analyses/Ryan/20200406'

photodat <- read.table(file.path(input_prefix,'Ryan/20200406/results_202004061414.txt'), header=T, sep='\t', quote = '', stringsAsFactors=FALSE)
photodat <- photodat[!is.na(photodat[,1]) & !is.na(photodat[,2]) & !is.na(photodat[,3]),]
photodat <- photodat[as.logical(photodat[,2]) & as.logical(photodat[,3]),]
rownames(photodat) <- paste0(rownames(photodat),'_',sapply(photodat[,1],function(x) strsplit(x,split='_',fixed=TRUE)[[1]][[4]]))

envvars <- c('bleached','parrotfish_scars','bivalve','Spirobranchus','Tridacna','boringurchin','other_polychaete','sponge','Halimeda','Turbinaria','Dictyota','Lobophora','cca','Galaxaura','Sargassum','unhealthy','turf_over','sediment','pigmentation','trematodiasis','gallcrabs','ascidians')

photonames <- gsub('-','',sub('^.*_OA000-','',rownames(photodat)))
photosamples <- unique(photonames)

photodatagg <- t(sapply(photosamples, function(x) {
    photosub <- photodat[photonames == x, c('porites_lumps','porites_ridges',envvars[envvars != 'bleached'])]
    if(sum(photonames == x) > 1) {
        apply(photosub,
              2,
              function(y) as.numeric(if(all(is.na(as.logical(y)))) NA else any(as.logical(y),na.rm=TRUE)))
    } else {
        return(as.numeric(as.logical(photosub)))
    }
}))

## host
allgenera <- unique(photodat$target)
ngen <- length(allgenera)

hostphotomat <- matrix(0, nrow = nrow(photodat), ncol = ngen)
colnames(hostphotomat) <- 1:ncol(hostphotomat)
rownames(hostphotomat) <- rownames(photodat)

hostphotomat[,1:ngen] <- sapply(allgenera, function(y) y == photodat$target)
colnames(hostphotomat)[1:ngen] <- paste('genus', allgenera, sep='_')


allmorphs <- sort(unique(paste(photodat$target,photodat$morphotype,sep='_'))[unique(photodat$morphotype) != ''])
nlevs <- length(allmorphs)
morphgroups <- sapply(allmorphs,function(x) strsplit(x,'_')[[1]][[1]])

hostphotomat <- cbind(hostphotomat, matrix(0,nrow=nrow(hostphotomat),ncol=nlevs))

hostphotomat[,ngen + (1:nlevs)] <- sapply(allmorphs, function(y) y == paste(photodat$target,photodat$morphotype,sep='_'))
colnames(hostphotomat)[ngen + (1:nlevs)] <- paste('morphotype', allmorphs, sep='_')

hostphotomat2 <- t(sapply(photosamples, function(x) if(sum(photonames == x) > 1) {apply(hostphotomat[photonames == x,],2,function(y) sum(y,na.rm=TRUE)/length(y))} else {hostphotomat[photonames == x,]}))

for(u in unique(morphgroups)) {
    relcols <- c(rep(TRUE,ngen), grepl(paste0('morphotype_',u),colnames(hostphotomat2)[(ngen+1):ncol(hostphotomat2)]))
    relrows <- hostphotomat2[,grepl(paste0('genus_',u),colnames(hostphotomat2))] == 1
    
    hostphotomat2[relrows,!relcols] <- NA
}

hostphotomat2 <- cbind(hostphotomat2, photodatagg[rownames(hostphotomat2),c('porites_lumps','porites_ridges')])

## env
bleachLevs <- unique(photodat[,'bleached'])[!is.na(unique(photodat[,'bleached'])) & unique(photodat[,'bleached']) != '']

envphotomat <- sapply(bleachLevs, function(y) y == photodat$bleached)
colnames(envphotomat) <- paste('bleached', bleachLevs, sep='_')
envphotomat2 <- t(sapply(photosamples, function(x) if(sum(photonames == x) > 1) {apply(envphotomat[photonames == x,],2,function(y) sum(y,na.rm=TRUE)/length(y))} else {envphotomat[photonames == x,]}))

envphotomat2 <- cbind(envphotomat2, photodatagg[,envvars[envvars != 'bleached']])
envphotomat2 <- envphotomat2[,apply(envphotomat2,2,sd,na.rm=TRUE) > 0]

photosOut <- cbind(hostphotomat2,envphotomat2)

write.table(photosOut, file = file.path(output_prefix, 'photo_annotations_aggregated_20201406.txt'), sep = '\t', quote=FALSE)
