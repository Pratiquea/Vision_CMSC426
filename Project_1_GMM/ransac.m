%{
Ransac flow:
1. Load the data and select s=2 number of points 
2. Fit the line using the two and then take the perp. distances
   from all points and reject those which are greater than 
    threshold distance 'd'
3. Meanwhile, keep logging the counts of the inliners.
4. Iterate through 1->3 , 'N' number of times. 
5. Get the max number of inliners in a particular iteration
%}
data=load('data/data1.mat');
ranSac(data);
function ranSac(data)
    %parameters
    s=2
    N=15
    d=16


    xlim=linspace(-200,150,1000);
    outliersx=0;
    outliersy=0;
    pts=data.pts;
    rawData=pts;
    rawDatax=transpose(pts(1,:));
    meanx=mean(rawDatax);
    rawDatay=transpose(pts(2,:));
    meany=mean(rawDatay);
    countArray=zeros(N,1);
    pointsArray=[0;0];
    inlinerStack=zeros(1,200);

    for i=1:N
        %select n points at random from the set
        sizeData=size(pts);
        randInd=randi([1 sizeData(:,2)],1,2);
        pt1=pts(:,randInd(1));
        pt2=pts(:,randInd(2));

        %getting line coff from these points
        coefficients = polyfit([pt1(1,:) pt2(1,:)],[pt1(2,:) pt2(2,:)], 1);
        a = coefficients (1);
        b = coefficients (2);
        cof=[a;-1];
        %getting the equation

        line=a*xlim+b;

        %gather all distances according to the equation ax-y+b=0
        dist2Line=abs(transpose(pts)*cof+b);

        %{check for the farthest point away
        [M,I]=max(dist2Line);
        pt=pts(:,I);

        %check for values that meet certain threshold d
        inliners=(dist2Line<d);
        size(inliners);
        [inlinCount,n]=size(dist2Line(inliners));
        inlinerStack=vertcat(inlinerStack,transpose(inliners));

        %store the counts 
        countArray(i)=inlinCount;
        pointsArray=horzcat(pointsArray,[a;b]);
    end

    %get max value of inlinecounts 
    [val,ind]=max(countArray);
    t=inlinerStack(ind+1,1:200);
    t=logical(t);
    outliersx=rawDatax(t);
    outliersy=rawDatay(t);

    %use the index to gather coff
    slope=pointsArray(1,ind+1);
    intercept=pointsArray(2,ind+1);
    fitLine=slope*xlim+intercept;
    %plot these two and visualize
    scatter(pts(1,:),pts(2,:),'k')
    hold on
    plot(outliersx,outliersy,'xg')
    %plot([pt1(1,:) pt2(1,:)],[pt1(2,:) pt2(2,:)],'r')
    plot(xlim,fitLine,'r')
    hold off
end

