########## Calculating feature ratios
#' @keywords internal
calculateRatios = function(currData){
  apply(currData,2,function(x){
    xRatios = matrix(x)%*%t(matrix(1/x))
    currIndexes = which(upper.tri(xRatios,diag = F),arr.ind = T)
    currRatios = xRatios[currIndexes]
    names(currRatios) = paste(names(x)[currIndexes[,1]], names(x)[currIndexes[,2]],sep="_")
    currRatios
  })
}

########## Building patient data
#' @keywords internal
dataCreation = function(trainData, sampleNames){
  lapply(unique(sampleNames), function(currSample){
    selectedIndexes = which(sampleNames == currSample)
    selectedTrajectoryBase = 1:length(selectedIndexes)
    trajectory= selectedTrajectoryBase/max(selectedTrajectoryBase)
    currData = trainData[,selectedIndexes]
    baseData = currData
    currDataNorm = t(apply(currData, 1, function(x){
      minVal = min(x)
      maxVal = max(x)
      if(max(x)==0){
        return(rep(0,length(x)))
      }else{
        (x - minVal)/(maxVal - minVal)
      }
    }))
    list(scaledData = currDataNorm, traj = trajectory, baseData = baseData, name = currSample, type = "Sample")
  })
}

########## Calculation of consensus trajectory
#' @keywords internal
createTrajectoryFromData = function(sampleData, sampleTraj, seed){
  pcaOfSample = stats::prcomp(t(sampleData[seed,]),scale. = T, center = T)
  distPCA = fields::rdist(pcaOfSample$x[,1:which(summary(pcaOfSample)$importance[3,]>0.9)[1]])
  sampleFixedTraj = cumsum(c(0,unlist(lapply(2:length(sampleTraj),function(i){
    abs(distPCA[i,i-1])
  }))))

  (sampleFixedTraj-min(sampleFixedTraj))/(max(sampleFixedTraj)-min(sampleFixedTraj))
}

########## Data smoothing
#' @keywords internal
computeNewData = function(sampleTraj, trajCond, dataToTransorm,winSz){
  dist2Others = outer(trajCond, -sampleTraj , "+")
  weightedData = exp(-(dist2Others^2)/(winSz^2))
  weightedData = weightedData/rowSums(weightedData)
  t(t(dataToTransorm) %*% t(weightedData))
}

########## Patient weigths for the alignment
#' @keywords internal
calculateSampleWeights = function(sample1Data,sample2Data){
  #nonNanFeatures = which(!is.nan(rowSums(sample1Data)) & !is.nan(rowSums(sample2Data)))
  #cors = stats::cor(sample1Data[nonNanFeatures,],sample2Data[nonNanFeatures,])
  cors = stats::cor(sample1Data,sample2Data)
  relevantCors = diag(cors)
  abs(mean(relevantCors[!is.na(relevantCors)]))
}

########## Pairwise alignment
#' @keywords internal
getAlignment = function(sample1, sample2){
  #sample1Data = t(sample1$scaledData)
  #sample2Data = t(sample2$scaledData)
  #selectedFeatures = which(!is.nan(colSums(sample1Data)) & !is.nan(colSums(sample2Data)))
  #disMatrix = fields::rdist(sample1Data[,selectedFeatures], sample2Data[,selectedFeatures])
  disMatrix = fields::rdist(t(sample1$scaledData), t(sample2$scaledData))
  dtw::dtw(disMatrix)
}

########## Model training (inner)
#' @keywords internal
multiAlign = function(listOfSamplesSmall, seed, numOfIter, no_cores){
  #### All pairwise alignments ####
  if(is.null(no_cores)){
    no_cores = max(1, parallel::detectCores() - 1)
  }
  message('Calculating all pairwise alignments:')
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("getAlignment", "listOfSamplesSmall", "seed"), envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = length(listOfSamplesSmall), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  `%dopar2%` <- foreach::`%dopar%`
  firstSampleInd = NULL
  allAlignments <- foreach::foreach(firstSampleInd = 1:length(listOfSamplesSmall), .options.snow = opts) %dopar2% {
    firstSample = listOfSamplesSmall[[firstSampleInd]]
    alignmentList = lapply(listOfSamplesSmall, function(secondSample){
      getAlignment(firstSample,secondSample)
    })
    setTxtProgressBar(pb, firstSampleInd)
    alignmentList
  }
  parallel::stopCluster(cl)
  close(pb)

  #### Consensus list ####
  message('Creating consensus list:')
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("calculateSampleWeights", "listOfSamplesSmall","getAlignment", "seed","createTrajectoryFromData","allAlignments"), envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = numOfIter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  `%dopar2%` <- foreach::`%dopar%`
  iterNum = NULL
  consensusList <- foreach::foreach(iterNum = 1:numOfIter, .options.snow = opts) %dopar2% {
    phylTree = phangorn::upgma(fields::rdist(sample(1:length(listOfSamplesSmall))))

    #### Couple joining ####
    nodesInTheTree = 1:(phylTree$Nnode+1)
    listOfSamplesForTree = listOfSamplesSmall
    for(i in 1:length(listOfSamplesForTree)){
      listOfSamplesForTree[[i]]$sampleScore = mean(sapply(1:length(allAlignments[[i]]),function(j){
        currAlignmentSteps = (allAlignments[[i]])[[j]]
        calculateSampleWeights(listOfSamplesSmall[[i]]$scaledData[,currAlignmentSteps$index1],
                               listOfSamplesSmall[[j]]$scaledData[,currAlignmentSteps$index2])
      }))
    }

    message('Creating a concensus sample...')
    edgeInformation = cbind(phylTree$edge,phylTree$edge.length)
    allFatherNodes = edgeInformation[,1]

    while(length(unique(nodesInTheTree))>1){
      currEdgeInformation = edgeInformation[edgeInformation[,1] %in% unique(allFatherNodes),]
      fatherOfKnownNodes = unique(allFatherNodes)[which(unlist(lapply(unique(allFatherNodes),function(i){
        length(which(currEdgeInformation[allFatherNodes==i,2] %in% nodesInTheTree))==2
      })))]
      knownNodes = unlist(lapply(fatherOfKnownNodes,function(i){currEdgeInformation[allFatherNodes==i,2]}))
      for(currFather in fatherOfKnownNodes){
        currNodes = currEdgeInformation[allFatherNodes==currFather,2]
        currSample1 = listOfSamplesForTree[[currNodes[1]]]
        currSample2 = listOfSamplesForTree[[currNodes[2]]]

        refSample = currSample1
        secSample = currSample2
        if(length(currSample1$traj)<length(currSample2$traj)){
          refSample = currSample2
          secSample = currSample1
          currNodes = rev(currNodes)
        }

        if(refSample$type=="Comb" | secSample$type=="Comb"){
          currAlignmentSteps = getAlignment(secSample,refSample)
        }else{
          currAlignmentSteps = (allAlignments[[currNodes[2]]])[[currNodes[1]]]
        }

        ## Weighting expression profiles
        refSamplePriority = refSample$sampleScore/(refSample$sampleScore+secSample$sampleScore)
        secSamplePriority = secSample$sampleScore/(refSample$sampleScore+secSample$sampleScore)

        newExp = refSample$scaledData[,currAlignmentSteps$index2]*refSamplePriority +
          secSample$scaledData[,currAlignmentSteps$index1]*secSamplePriority
        newTraj = seq(0,1,length.out = dim(newExp)[2])

        newBase = refSample$baseData[,currAlignmentSteps$index2]*refSamplePriority +
          secSample$baseData[,currAlignmentSteps$index1]*secSamplePriority

        fixedTraj = newTraj
        alignQuality = calculateSampleWeights(refSample$scaledData[,currAlignmentSteps$index2],
                                              secSample$scaledData[,currAlignmentSteps$index1])


        listOfSamplesForTree[[currFather]] = list(scaledData = newExp, traj = fixedTraj, name = currFather, type = "Comb", sampleScore = alignQuality, baseData = newBase)
      }
      nodesInTheTree = c(nodesInTheTree[!(nodesInTheTree %in% knownNodes)],fatherOfKnownNodes)
      allFatherNodes = allFatherNodes[!(allFatherNodes %in% fatherOfKnownNodes)]
    }

    #### Interpreting results ####
    fullAlignedSample = listOfSamplesForTree[[nodesInTheTree]]
    fullAlignedSample$traj = createTrajectoryFromData(fullAlignedSample$scaledData,fullAlignedSample$traj, seed)

    setTxtProgressBar(pb, iterNum)
    fullAlignedSample
  }
  parallel::stopCluster(cl)
  close(pb)
  consensusList
}

#' Training the TimeAx model
#'
#' This function initiate model training using the TimeAx algorithm - performing a multiple trajectory alignment (MTA) on time-series datasets of individuals each of which is considered as an individual partial trajectory.
#'
#' @param trainData A matrix containing profiles (columns) of omics measurments (rows) from multiple individuals and different time points. For omics data it is better to use raw values instead of normalized ones. Profiles for each individual should be ordered by chronological time.
#' @param sampleNames A vector containing the individual identity of each sample in the train data.
#' @param ratio Boolean parameter determining whether the model will be based on the ratios between feature values or the raw values. The default is to use ratios (TRUE) as this allows the model to overcome technical batch effects. For data types with many zeros or discrete values, the raw data should be used (selecting FALSE).
#' @param numOfIter Number of consensus trajectories. The default is 100.
#' @param numOfTopFeatures Length of the conserved-dynamics-seed of features. If trainData has less features, all of them will be selected as seed. The default is 50.
#' @param seed The conserved-dynamics-seed. If provided, the alignment process will be conducted based on these features. The default is NULL.
#' @param no_cores A number for the amount of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.
#' @return A TimeAx model consists of:
#' \item{consensusList}{List of consensus trajectories.}
#' \item{seed}{The conserved-dynamics-seed that was used in the alignment process.}
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Training the model
#' model = modelCreation(DataUBC,UBCSamples, no_cores = 2)
#'
#' \dontrun{
#'
#' }
#' @export
#' @importFrom "utils" "setTxtProgressBar"
#' @importFrom "stats" "sd" "var"
#' @importFrom "grDevices" "chull"
modelCreation = function(trainData, sampleNames, ratio = T, numOfIter = 100, numOfTopFeatures = 50 ,seed = NULL, no_cores = NULL){
  listOfSamples = dataCreation(trainData, sampleNames)

  # Checking that features are variable across all patients
  # if(length(which(is.nan(do.call(rbind,lapply(listOfSamples,function(x){rowSums(x$scaledData)})))))>0){
  #   message('Some features do not vary across patients time points!')
  #   return(NULL)
  # }
  if(max(apply(do.call(rbind,lapply(listOfSamples, function(X){
    apply(X$scaledData,1,function(x){
      length(unique(x))
    })
  })),2,min))==1){
    message('None of the features vary across patients time points!')
    return(NULL)
  }

  # Checking if there are negative values
  if(min(trainData)<0){
    message('Negative values, switching ratios off!')
    ratio = F
  }

  if(is.null(seed)){
    if(dim(trainData)[1]<numOfTopFeatures){
      message(paste('Seed size switched to ', dim(trainData)[1],sep = ""))
      seed = row.names(trainData)
    }else{
      seed = detectSeed(listOfSamples, sampleNames, numOfTopFeatures = numOfTopFeatures, no_cores = no_cores)
    }
  }


  listOfSamplesSmall = lapply(listOfSamples,function(currSample){
    currSampleNew = currSample
    currSampleNew$scaledData = currSample$scaledData[seed,]
    currSampleNew$baseData = currSample$baseData[seed,]
    currSampleNew
  })

  consensusList = multiAlign(listOfSamplesSmall, seed, numOfIter = numOfIter, no_cores = no_cores)

  message('Creating a TimeAx model')
  consensusList = lapply(consensusList,function(x){
    if(ratio){
      currData = calculateRatios(x$baseData)
    }else{
      currData = x$baseData
    }
    list(baseData = currData, traj = x$traj)
  })

  #### Output ####
  message('Model created')
  list(consensusList = consensusList, seed = seed, ratio = ratio)
}

#' Selecting conserved-dynamics-seed features
#'
#' @param trainData A matrix containing profiles (columns) of omics measurments (rows) from multiple individuals and different time points. For omics data it is better to use raw values instead of normalized ones. Profiles for each individual should be ordered by chronological time.
#' @param sampleNames A vector containing the individual identity of each sample in the train data.
#' @param numOfTopFeatures Length of the conserved-dynamics-seed of features. If trainData has less features, all of them will be selected as seed. The default is 50.
#' @param topGenes Number of initial high variable features to be considered for the seed selection. The default is 4000.
#' @param numOfIterations Number of different random sample selections for the calculation. The default is 20.
#' @param percOfSamples Fraction of samples from each individual, selected in each sample selection. The default is 0.8
#' @param no_cores A number for the amount of cores which will be used for the analysis. The default (NULL) is total number of cores minus 1.
#' @return A list including:
#' The conserved-dynamics-seed. A list of feature, suitable for the model training, ordered from best to worst.
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Selecting conserved-dynamics-seed features
#' seed = detectSeed(DataUBC,UBCSamples, no_cores = 2)
#' @export
detectSeed = function(trainData, sampleNames, numOfTopFeatures = 50, topGenes = 4000, numOfIterations = 20, percOfSamples = 0.8, no_cores = NULL){
  listOfSamples = trainData
  if(!is.null(dim(trainData))){
    if(dim(trainData)[1]<numOfTopFeatures){
      message(paste('Seed size switched to ', dim(trainData)[1],sep = ""))
      return(row.names(trainData))
    }
    listOfSamples = dataCreation(trainData, sampleNames)
  }

  baseDataList = lapply(listOfSamples,function(x){x$scaledData})
  minSize = round(min(unlist(lapply(baseDataList,function(x){dim(x)[2]})))*percOfSamples)

  message('Initial feature selection')
  mutualGeneNames = row.names(baseDataList[[1]])

  # Removing features with only one unique value in some of the patients #
  mutualGeneNames = mutualGeneNames[apply(do.call(rbind,lapply(baseDataList, function(X){
    apply(X,1,function(x){
      length(unique(x))
    })
  })),2,min)>1]

  if(length(mutualGeneNames)<(2*topGenes)){topGenes = floor(length(mutualGeneNames)/2)}

  # Combining the data of all subjects into one big data frame #
  baseDataList = lapply(baseDataList,function(x){x[mutualGeneNames,]})
  dataForHighVar = as.data.frame(t(do.call(cbind,baseDataList)))

  # Filtering features with many non-unique values #
  maxUniques = dim(dataForHighVar)[1]+2-2*length(baseDataList) # Counting for all subjects together and removing the first and last samples which are always 0 and 1
  minAllowedUniques = 0.8*maxUniques
  genesUniques = lengths(lapply(dataForHighVar, unique))
  selectedGeneIndexes = which(genesUniques>minAllowedUniques)
  if(length(selectedGeneIndexes)<2*topGenes){
    selectedGeneIndexes = order(genesUniques,decreasing = T)[1:(2*topGenes)]
  }

  # Filtering features with low standard deviation #
  genesSD = apply(dataForHighVar[,selectedGeneIndexes],2,sd)
  selectedGeneIndexes = selectedGeneIndexes[order(genesSD,decreasing = T)[1:topGenes]]
  mutualGeneNames = mutualGeneNames[selectedGeneIndexes]

  # Filtered subject list #
  baseDataList = lapply(baseDataList,function(x){x[mutualGeneNames,]})

  message('Choosing conserved-dynamics-seed:')
  if(is.null(no_cores)){
    no_cores = min(numOfIterations,max(1, parallel::detectCores() - 1))
  }
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("mutualGeneNames", "minSize","baseDataList"), envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = numOfIterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  `%dopar2%` <- foreach::`%dopar%`
  iteration = NULL
  geneScoreList = foreach::foreach(iteration = 1:numOfIterations, .options.snow = opts) %dopar2% {
    sampledData = t(do.call(rbind,lapply(baseDataList, function(currSample){
      selectedCells = sort(sample(1:dim(currSample)[2],minSize))
      newData = currSample[,selectedCells]
      colnames(newData) = 1:dim(newData)[2]
      newData
    })))

    setTxtProgressBar(pb, iteration)

    sapply(1:length(mutualGeneNames), function(i){
      currMatrix = sampledData[,seq(i,dim(sampledData)[2],length(mutualGeneNames))]
      mean(stats::cor(currMatrix,method = "spearman"))
    })
  }
  parallel::stopCluster(cl)
  close(pb)

  geneScore = colMeans(do.call(rbind,geneScoreList))
  mutualGeneNames[order(geneScore,decreasing = T)[1:numOfTopFeatures]]
}

#' Infer pseudotime for new samples bassed on the TimeAx model
#'
#' @param model A TimeAx model.
#' @param testData A matrix containing profiles (columns) of features measurments (rows). Data should provided in similar scales (preferably, non-normalized) as the train data. Seed genes that are missing in the test data will be excluded from the prediction.
#' @param no_cores A number for the amount of cores which will be used for the analysis. The default (NULL) is total number of cores minus 1.
#' @param seed The conserved-dynamics-seed. If provided, the prediction process will be conducted based on these features. Use the model's seed by keeping the the default value of NULL.
#' @param sampleNames Used for the robustness analysis. Always keep as NULL.
#' @return A prediction list consists of:
#' \item{predictions}{The final pseudotime position for each sample.}
#' \item{uncertainty}{An uncertainty score for each position. Lower scores means higher certainty.}
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Training the model
#' model = modelCreation(DataUBC,UBCSamples, no_cores = 2)
#'
#' # Inferring pseudotime positions
#' pseudotimeStats = predictByConsensus(model,DataUBC, no_cores = 2)
#' pseudotime = pseudotimeStats$predictions
#' uncertainty = pseudotimeStats$uncertainty
#' @export
predictByConsensus = function(model, testData, no_cores = NULL, seed = NULL, sampleNames = NULL){
  if(is.null(seed)){
    seed = intersect(model$seed, row.names(testData))
  }

  if(is.null(no_cores)){
    no_cores = max(1, parallel::detectCores() - 1)
  }

  ratio = model$ratio

  testData = testData[seed,,drop = F]
  if(ratio){
    testData = calculateRatios(testData)
  }

  cleanModel = lapply(model$consensusList,function(x){
    #currData = calculateRatios(x$baseData[seed,])
    currData = x$baseData[row.names(testData),]
    list(data = currData, traj = x$traj)
  })

  if(is.null(sampleNames)){
    message('Predicting samples pseudotime positions:')
  }else{
    message('Predicting robustness pseudotime positions:')
  }
  cl<-parallel::makeCluster(no_cores)
  parallel::clusterExport(cl=cl, varlist=c("seed", "testData","sampleNames","computeNewData","cleanModel"), envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(min = 1, max = length(cleanModel), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  `%dopar2%` <- foreach::`%dopar%`
  consensusInd = NULL
  predictionStats <- foreach::foreach(consensusInd = 1:length(cleanModel), .options.snow = opts) %dopar2% {
    sampleConsensus = cleanModel[[consensusInd]]
    refNorm = sampleConsensus$data
    testDataNorm = testData

    pairsToUse = which(rowSums(testDataNorm)!=Inf & !is.nan(rowSums(testDataNorm)) &
            rowSums(refNorm)!=Inf & !is.nan(rowSums(refNorm)))

    refNorm = refNorm[pairsToUse,]
    testDataNorm = testDataNorm[pairsToUse,,drop=F]

    refMeans = rowMeans(refNorm)
    refTestRatio = stats::median(refMeans/rowMeans(testDataNorm))
    refSD = apply(refNorm,1,sd)
    refNorm = (refNorm - refMeans)/refSD

    testDataNorm = (testDataNorm*refTestRatio-refMeans)/refSD

    selectedPairs = row.names(testDataNorm)
    outDist = rowMeans(testDataNorm)
    outDistP = 2*stats::pnorm(abs(outDist), lower.tail = F)
    selectedPairs = names(which(outDistP>0.1))
    if(length(selectedPairs)<length(seed)){
      selectedPairs = row.names(testDataNorm)
    }

    refNorm = refNorm[selectedPairs,]
    testDataNorm = testDataNorm[selectedPairs,,drop=F]

    if(is.null(sampleNames)){
      corMatrix = stats::cor(refNorm,testDataNorm,method = "spearman")
      corMatrix = computeNewData(sampleConsensus$traj,sampleConsensus$traj,corMatrix,0.1)
      currMaxIndexes = apply(corMatrix,2,which.max)
      prediction = sampleConsensus$traj[currMaxIndexes]
    }else{
      prediction = rep(NA,length(sampleNames))
      for(currSample in unique(sampleNames)){
        currIndexes = which(currSample == sampleNames)
        testSample = testDataNorm[,currIndexes,drop=F]
        corMatrix = t(stats::cor(as.matrix(testSample),refNorm,method = "spearman"))
        corMatrix = 1 - computeNewData(sampleConsensus$traj,sampleConsensus$traj,corMatrix,0.1)
        startMatrix = cbind(0,t(as.matrix(corMatrix[,1])))
        matrixList = NULL
        if(dim(corMatrix)[2]>1){
          matrixList = lapply(2:dim(corMatrix)[2], function(i){
            currConsMatrix = t(matrix(corMatrix[,i],nrow = dim(corMatrix)[1],ncol = dim(corMatrix)[1]))
            currConsMatrix[lower.tri(currConsMatrix)] <- 0
            currConsMatrix
          })
        }
        endMatrix = rbind(matrix(1, ncol = 1, nrow = dim(corMatrix)[1]),0)

        namesForBigMatrix = c("Start", paste("Cons",1:dim(corMatrix)[1],1,sep="_"))
        if(dim(corMatrix)[2]>1){
          namesForBigMatrix = c(namesForBigMatrix, unlist(lapply(2:dim(corMatrix)[2], function(i){
            paste("Cons",1:dim(corMatrix)[1],i,sep="_")
          })))
        }
        namesForBigMatrix = c(namesForBigMatrix,"End")

        matrixList = c(list(startMatrix),matrixList,list(endMatrix))
        bigMatrix = Matrix::.bdiag(matrixList)
        colnames(bigMatrix) = row.names(bigMatrix) = namesForBigMatrix

        dfForGraph = data.frame(from = namesForBigMatrix[bigMatrix@i+1], to = namesForBigMatrix[bigMatrix@j+1], cost = bigMatrix@x)
        modelGraph = cppRouting::makegraph(dfForGraph,directed = T)
        shortestPath<-cppRouting::get_path_pair(modelGraph,from="Start",to="End")
        pathNodes = rev(shortestPath$Start_End)
        finalIndexesForPrediction = as.numeric(sapply(pathNodes[2:(length(pathNodes)-1)],function(x){unlist(strsplit(x,"_"))[2]}))
        currPrediction = sampleConsensus$traj[finalIndexesForPrediction]

        setTxtProgressBar(pb, consensusInd)
        prediction[currIndexes] = currPrediction
      }
    }
    prediction
  }
  parallel::stopCluster(cl)
  close(pb)

  predictionMatrix = do.call(rbind,predictionStats)
  finalPredictions = colMeans(predictionMatrix)
  sampleUnCertainty = apply(predictionMatrix,2,sd)
  list(predictions = finalPredictions, uncertainty = sampleUnCertainty)
}

#' Calculate a robustness score for the TimeAx model
#'
#' @param model A TimeAx model.
#' @param trainData The matrix containing profiles (columns) of omics measurments (rows), which was used to train the model.
#' @param sampleNames A vector containing the individual identity of each sample in the train data. Same vector as used in the training.
#' @param pseudo The output list of predictByConsensus. If not provided (NULL), pseudotime will be inferred by this function.
#' @param no_cores A number for the amount of cores which will be used for the analysis. The default (NULL) is total number of cores minus 1.
#' @return A robustness list consists of:
#' \item{robustnessPseudo}{Robustness pseudotime positions for all samples.}
#' \item{score}{TimeAx robustness score for the model.}
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Training the model
#' model = modelCreation(DataUBC,UBCSamples,no_cores = 2)
#'
#' # Inferring pseudotime positions
#' robustnessStats = robustness(model,DataUBC,UBCSamples,no_cores = 2)
#' robustnessPseudo = robustnessStats$robustnessPseudo
#' robustnessScore = robustnessStats$score
#' @export
robustness = function(model, trainData, sampleNames, pseudo = NULL, no_cores = NULL){
  if(is.null(pseudo)){
    pseudo = predictByConsensus(model, trainData, no_cores = no_cores)$predictions
  }else{
    pseudo = pseudo$predictions
  }
  pseudoRobust = predictByConsensus(model, trainData, no_cores = no_cores, sampleNames = sampleNames)$predictions
  list(robustnessPseudo = pseudoRobust, score = stats::cor(pseudoRobust,pseudo))
}

#' Calculate K-fold cross valudation for the TimeAx model
#'
#' @param model A TimeAx model.
#' @param trainData The matrix containing profiles (columns) of omics measurments (rows), which was used to train the model.
#' @param sampleNames A vector containing the individual identity of each sample in the train data. Same vector as used in the training.
#' @param k Number of folds to use.
#' @param no_cores A number for the amount of cores which will be used for the analysis. The default (NULL) is total number of cores minus 1.
#' @return A K-fold cross validation score.
#' @references
#' Submitted
#' @examples
#' data(UBCData)
#'
#' # Training the model
#' model = modelCreation(DataUBC,UBCSamples,no_cores = 2)
#'
#' # Calculating K-fold cross validation score
#' kFoldScore = kFold(model,DataUBC,UBCSamples, k = 5)
#' @export
kFold = function(model, trainData, sampleNames, k, no_cores = NULL){
  message("Calculate global disease pseudotime")
  pseudoAll = predictByConsensus(model,trainData,no_cores = no_cores)$predictions

  message("\nCalculate K-fold pseudotime")
  patientSplit = split(unique(sampleNames), ceiling(seq_along(unique(sampleNames)) / ceiling(length(unique(sampleNames))/k)))
  pseuedoK = do.call(rbind, lapply(1:k, function(i){
    message(paste('\nFold',i,sep = " "))
    testPatients = patientSplit[[i]]
    testIndexes = which(sampleNames %in% testPatients)
    trainPatients = unique(sampleNames)[!(unique(sampleNames) %in% testPatients)]
    currModel = modelCreation(trainData[,sampleNames %in% trainPatients],sampleNames[sampleNames %in% trainPatients],seed = model$seed,ratio = model$ratio, no_cores = no_cores)
    cbind(testIndexes,predictByConsensus(currModel,trainData[,sampleNames %in% testPatients], no_cores  = no_cores)$predictions)
  }))
  pseuedoKSorted = rep(1,dim(pseuedoK)[2])
  pseuedoKSorted[pseuedoK[,1]] = pseuedoK[,2]
  list(Score = cor(pseudoAll,pseuedoKSorted), kFoldPseudo = pseuedoKSorted)
}


#' UBC RNA-seq data.
#'
#' @format A matrix with 27 profiles (columns) of 5000 genes.
"DataUBC"

#' UBC sample labels.
#'
#' @format A vector with 27 sample labels for 5 individuals.
"UBCSamples"
