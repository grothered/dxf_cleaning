##
dta_in=read.csv('tmp10.txt',header=F)

dta=matrix(NA,nrow=dim(dta_in)[1],ncol=3)
dta[,1]=as.numeric(as.character(dta_in[,1]))
dta[,2]=as.numeric(as.character(dta_in[,2]))
dta[,3]=as.numeric(as.character(dta_in[,3]))

# Check for valid input data
if(sum(is.na(dta[,1]))>0){
    print('ERROR: Non numeric values in the input data')
    print('These are:')
    tmp=which(is.na(dta[,1]))
    print(cbind(tmp, dta_in[tmp,]))
    stop()
    print('You need to go into the files and manually fix this')
}


ycopy=dta[,3]*NA
xcopy=dta[,2]*NA
for(i in 1:length(dta[,3])){
    ycopy[i] = sum(abs(dta[,3]-dta[i,3])<2.0e-00)
    xcopy[i] = sum(abs(dta[,2]-dta[i,2])<2.0e-00)
}
#
## Find points which might be part of xy groups They will have the
## corresponding x/y point directly above/below them
v=which((xcopy>1)&(ycopy>3))
#
#plot( dta[v,2:3],cex=0.2 )

# Idea -- Find pairs of points with the same x value, check if they both have the same number of y neighbours, and see if this forms a cluster.

count=0
yz=dta[,3]
xz=dta[,2]

# Parameters
yrange=400
ytol=5
xtol=2
taken=dta[,3]*0 # This will record points that we have used

for(i in v){
    #print(c('Next v', i))
    # Find points with similar y values
    row1=which(abs(yz-dta[i,3])<ytol)
    if(length(row1)==0) next
    # Exclude points which have already been taken by other cross-sections
    row1=row1[which(taken[row1]==0)]
    if(length(row1)==0) next
    
    # Sort them according to increasing x values
    row1_order=sort(dta[row1,2],index.return=T)$ix
    row1=row1[row1_order]

    # Make sure they are all within 400 x-units of each other
    row1_xvals=dta[row1,2]
    row1_xdist=abs(row1_xvals - row1_xvals[1])
    row1=row1[which(row1_xdist<400)]
    #if(length(row1)<=3){
    #    print(c('Too few similar y values in row1', length(row1)))
    #    next 
    #}
    row1_xvals=dta[row1,2]


    # Find points with nearly the same x values as points in row1
    row2=c()
    #print("Making row2 ...")
    for(j in row1_xvals){
        row2=c(row2, which((xz>j-xtol)&(xz<j+xtol)))
    }

    # Exclude points already in row1 from row2
    row2=setdiff(row2,row1)
    #if(length(row2)<length(row1)){
    #     print(c('Too few row2 values, exclusion step 1', length(row2), length(row1)))
    #     next
    #}

    
    # Sort the row2 values
    row2_order=sort(dta[row2,2],index.return=T)$ix
    row2=row2[row2_order]
    row2_xvals=dta[row2,2]
    row2_yvals=dta[row2,3]
    row1_yvals=dta[row1,3]

    # Figure out which row2 values correspond to row1 values
    row2_to_row1=row2*NA
    for(j in 1:length(row2)){
        row2_to_row1[j] = which.min(abs(row2_xvals[j] - row1_xvals))
    }

    # Now look for multiple row2 points matching a single row1 point
    repeated_values = row2_to_row1[which(duplicated(row2_to_row1))]
    repeated_indices=which(row2_to_row1 %in% repeated_values)
    non_repeated_indices=setdiff(1:length(row2), repeated_indices) #setdiff(row2_to_row1, row2_to_row1[which(repeated_points)])
    #print(c('Repeated points ', repeated_points))
    #print(c('Non-repeated_points', non_repeated_points))

    if(length(non_repeated_indices)>0){
        typical_y_dist=median( (dta[row2,3]-dta[row1[row2_to_row1],3] )[non_repeated_indices] )
    }else{
        typical_y_dist=0
    }
    #print(c('Typical y dist is', typical_y_dist))

    if(length(repeated_values)>0){
        # Scroll through the repeated points, and remove all but the nearest neighbour
        exclude=c()
        for(j in repeated_values){
            f=which(row2_to_row1==j)
            xtmp=row1_xvals[j]
            ytmp=row1_yvals[j]
            
            # Keep the one which is the closest to the typical distance from the row1 point
            keep_test = which.min( abs( (row2_yvals[f]-ytmp) -typical_y_dist) )
            exclude = c(exclude, f[-keep_test])
        }
        # Cut the exclusion points out of row 2
        row2=row2[-exclude]
        row2_xvals=dta[row2,2]
        row2_yvals=dta[row2,3]
    }

    if(length(row2)<length(row1)){
        row1_keep=union(row2_to_row1,row2_to_row1)
        row1=row1[row1_keep]

    }

    if(length(row1)==0) next

    # Error check
    if(length(row2)!=length(row1)){
        print("ERROR - length row2 is not equal to length row1")
        print(row1)
        print(row2)
        stop()
    }

    # Now remove points which clearly still shouldn't be here
    typical_y_dist=median(dta[row1,3]-dta[row2,3])
    exclude=c()
    for(i in 1:length(row2)){
        y_row1=dta[row1[i],3]
        y_row2=dta[row2[i],3]
        ydist=y_row1-y_row2
        if(sign(ydist)!=sign(typical_y_dist)){
            exclude=c(exclude,i)
            next
        }
        if(abs(ydist)>1.5*abs(typical_y_dist)){
            exclude=c(exclude,i)
            next
        }
        if(abs(ydist)<1/1.5*abs(typical_y_dist)){
            exclude=c(exclude,i)
            next
        }

    }
    if(length(exclude>0)){
        row1=row1[-exclude]
        row2=row2[-exclude]
    } 

    if(length(row1)<6) next

    ### PREPARE DATA FOR OUTPUT

    # Find which row is top and which row is bottom 
    row1_y=median(dta[row1,3])
    row2_y=median(dta[row2,3])
    if(row2_y>row1_y){
        if(min(dta[row2,3])<max(dta[row1,3])){
            print('ERROR -- y values are overlapping')
            print(cbind(dta[row2,2:3], dta[row1,2:3]))
            stop()
        }
        xvalues=dta[row1,1]
        yvalues=dta[row2,1]

    }else{
        if(min(dta[row2,3])>max(dta[row1,3])){
            print('ERROR -- y values are overlapping')
            print(cbind(dta[row2,2:3], dta[row1,2:3]))
            stop()
        }
        xvalues=dta[row2,1]
        yvalues=dta[row1,1]

    }       

    ## X values on top 
    section=cbind(xvalues,yvalues)
    if(max(abs(section[,1]-sort(section[,1])))!=0){
        print('ERROR: section values are not increasing')
        print(section)
        stop()
    }else{
        print('SUCCESS? Section values are:')
        print(section)
        dev.new()
        plot(section[,1],section[,2],t='o',asp=4)
    }

    # Cut these points out of the list we search through
    taken[row1]=1
    taken[row2]=1
   
}



