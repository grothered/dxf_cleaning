## CODE TO EXTRACT CROSS_SECTIONS FROM THE TEXT LABELS IN THE DFX FILES

# Pasig River Cross-sections
all_files=dir(pattern="PR*.*.txt")

# Parameters
yrange=400 # The code looks for neighbouring points within this y-range to form the cross-section
ytol=5 # The code considers points with y-values that differ by ytol to be 'nearly horizontally aligned'
xtol=2 # The code considers points with x-values that differ by xtol to be 'nearly vertically aligned'


manually_clean<-function(dta_in,a_file){
    # This function allows the user to manually place edits into the data.  The
    # advantage of doing it this way is that the original data is cleanly
    # edited
    if(a_file=="PR00150-00250.txt"){
        # Correct data entry error
        ind=which((dta_in[,1]==96.17)) #&(dta_in[,2]== 86.476591)&(dta_in[,3]== -355.6864))
        print(c('Replacing ', ind, ' in ', a_file) )
        dta_in[ind,1]="86.87"
    }

    if(a_file=="PR00745-00800.txt"){
        # This section has a bridge, and the bridge elevations are mixed with the cross-section ones
        print(c('REMOVING SECTION WITH BRIDGE FROM ', a_file, '. You will have to treat this one manually '))
        dta_in = dta_in[dta_in[,3]<500,]
    }

    if(a_file=="PR01750-01850.txt"){
        # In this one, my script seems to have incorrectly read the output,
        # because there is a space after the minus sign in the character (unusual).
        ind=which(as.character(dta_in[,1])=="- ")
        print(c('Replacing ', ind, 'in ', a_file))        
        #stop()
        dta_in[ind,1]="-5.73"
    }


    if(a_file=="PR02350-02450.txt"){
        ind = which( (dta_in[,1]==13.80)|(dta_in[,1]==16.35))
        print(c('Removing ', ind , ' in ', a_file))
        dta_in=dta_in[-ind,]

        ind = which(dta_in[,1]==-43.95)
        print(c('Correcting ', ind , ' in ', a_file))
        dta_in[ind,1]=-33.95
    }    

    if(a_file=="PR03092-03200.txt"){
        print(c('REMOVING SECTION WITH BRIDGE FROM ', a_file, '. You will have to treat this one manually '))
        ind=which(dta_in[,3]>400)
        dta_in=dta_in[-ind,]
    }
        
    if(a_file=="PR03850-03950.txt"){
        ind1=which(dta_in[,1]==-67.82)
        ind2=which(dta_in[,1]==-67.78)
        print(c('Re-ordering points representing thin wall in ', a_file))
        #tmp=dta_in[ind1,]
        dta_in[ind1,1]=-67.78 #dta_in[ind2,1]
        dta_in[ind2,]=-67.82 #tmp

        ind=which(dta_in[,1]==63.18)
        ind2=ind[which.max(dta_in[ind,2])]
        print(c('Correcting ', ind2, ' in ', a_file))
        dta_in[ind2,1]=70.78
   
        ind=which(dta_in[,1]==64.78)
        print(c('Correcting ', ind, ' in ', a_file))
        dta_in[ind,1]=74.78 

    }


    if(a_file=="PR04000-04100.txt"){
        ind = which(dta_in[,1]==66.59)
        print(c('Correcting ', ind, ' in ', a_file))
        dta_in[ind,1]=59.69
    }

    if(a_file=="PR04600-04700.txt"){
        ind=which(dta_in[,1]==69.92)
        ind2=ind[which.max(dta_in[ind,2])]
        print(c('Correcting ', ind2, ' in ', a_file))
        dta_in[ind2,1]=70.20
    }


    # Return the modified data set
    return(dta_in)
}


# Loop through all files
for(a_file in all_files){

    # Read in the data
    dta_in=read.csv(a_file,header=F)
    dta_in[,1]=as.character(dta_in[,1])

    dta_in=manually_clean(dta_in,a_file)

    # Coerce the data into a suitable format
    dta=matrix(NA,nrow=dim(dta_in)[1],ncol=3)
    dta[,1]=as.numeric(as.character(dta_in[,1]))
    dta[,2]=as.numeric(as.character(dta_in[,2]))
    dta[,3]=as.numeric(as.character(dta_in[,3]))

    # Check for non-numeric characters in the input data
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
        ycopy[i] = sum(abs(dta[,3]-dta[i,3])<ytol)
        xcopy[i] = sum(abs(dta[,2]-dta[i,2])<xtol)
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
            print(c('ERROR: section values are not increasing in', a_file, 'around y= ', median(dta[row1,3])))
            print(section)
            stop()
        }else{
            print('SUCCESS? Section values are:')
            print(section)
            #dev.new()
            plot(section[,1],section[,2],t='o',asp=4)
        }

        # Cut these points out of the list we search through
        taken[row1]=1
        taken[row2]=1
   
    }

}



