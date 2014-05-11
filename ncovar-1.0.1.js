////////////////////////////////////////////////////////////////////////////////////////////////
//
// ncovar-1.0.1.js
//
// JavaScript function to calculate the test statistic and p-values of the NCoVaR test from Diks and Wolski (2014)
// 
// main function: ncovar(dataA,dataB,VaR,safe,bandwidth,forward,backward,scenario,variance,kernel,norm) { }
// dataA: time series for variable A
// dataB: time series for variable B (has to be the same length)
// VaR: Value-at-Risk level (for instance 0.05)
// safe: safe quantile level (for instance 0.5)
// bandwidth: bandwidth value
// forw: end of the forward-looking variable (0 for contemporary)
// backward: number of lags (0 for contemporary)
// scenario: "sc" for scenario 1 (tail to tail) and "sc2" for scenario 2 (median to tail)
// variance: variance estimator ("standard")
// kernel: kernel function, default is the "rectangular"
// norm: data normalisation ("standard" for standard formal and "uniform" for uniform transformation)
//
////////////////////////////////////////////////////////////////////////////////////////////////

//main function
function ncovar(dataA,dataB,VaR,safe,bandwidth,forward,backward,scenario,variance,kernel,norm) {
//data normalisation
if(norm=="standard") {
	var A=normalise(dataA);
	var B=normalise(dataB);
}
else {
	var A=uniform(dataA);
	var B=uniform(dataB);
}

//get quantiles depending on scenario considered
var qrA=getQuantile(VaR,A); //qxiq
var qsA=getQuantile(safe,A); //qxis
var qrB=getQuantile(VaR,B); //qrz
if(scenario=="sc1") {
	var qsB=qrB;
}
else if(scenario=="sc2") {
	var qsB=getQuantile(safe,B); //
}

//calculate T-values
return(tVal(A,B,forwStart,forward,backward,qrB,qsB,qrA,qsA,bandwidth,variance,kernel));
}

//calculate test statistic
function tVal(dataA,dataB,forwStart,forw,back,qz,qxj,qxir,qxis,eps,variance,kernel) {
var n=dataA.length;
var K1=[];

for(var k=back; k<(n-forw); k++) {
	var CIzr=0, CIzs=0, CIr=0, CIs=0;
	var Izr=0, Izs=0, Ir=0, Is=0;
	
	//calculate kernels (I) for index k
	
	switch(kernel) {
		default: 
			Ir=normKernel(qxj,dataB.slice(k-back,k+1),eps)*normKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*normKernel(qz,dataB.slice(k,k+forw+1),eps);
			Is=normKernel(qxj,dataB.slice(k-back,k+1),eps)*normKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*normKernel(qz,dataB.slice(k,k+forw+1),eps);
	}
	
	var m=back;
	while (m<=k) {
		//calculate correlation kernels (CI) for index m
		switch(kernel) {
			default:
				CIr=CIr+normKernel(qxj,dataB.slice(m-back,m+1),eps)*normKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+normKernel(qxj,dataB.slice(m-back,m+1),eps)*normKernel(qxir,dataA.slice(m-back,m+1),eps)*normKernel(qz,dataB.slice(m,m+forw+1),eps);
				CIs=CIs+normKernel(qxj,dataB.slice(m-back,m+1),eps)*normKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+normKernel(qxj,dataB.slice(m-back,m+1),eps)*normKernel(qxis,dataA.slice(m-back,m+1),eps)*normKernel(qz,dataB.slice(m,m+forw+1),eps);
		
		}
		m++;
		
	}
	//Correlation Integral
	K1.push(Izr*CIs-Ir*CIzs+CIzr*Is-CIr*Izs);
	
}
var forL=forw-forwStart+1;
var T=Math.pow(2*eps,(-forL-2*back-2*back))/Math.pow(n,2)*sum(K1);
var kernel=[];
for(var i=0; i<K1.length; i++) {
	kernel.push((Math.pow(2*eps,(-forL-2*back-2*back)))/n*(K1[i]));
}

//variance of the test statistics
var tVar=myVar(kernel);

var T_VAL=T*Math.sqrt(n)/Math.sqrt(4*tVar);
return T_VAL;
}

//rectangular kernel
function normKernel(point,vector,eps) {
var value=1;
	for(var i=0; i<vector.length; i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(p<1) {
			control=1;
		}
		value=Number(value*control);
	}
return value;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// data normalisation functions
////////////////////////////////////////////////////////////////////////////////////////////////

//data normalisation (standard normal)
function normalise(data) {
	if (data.length == 0) return null;
	var mean_value = mean(data); 
	var sd_value = sd(data);
	var tempData = [];
	for(var i=0; i<data.length; i++) {
		tempData.push((data[i]-mean_value)/sd_value);
	}
	return tempData;
}

//insertion sort for uniform transformation
function insertionSort(X) {
	var M=X.length;
	var I=[], S=[];
	var R=0, r=0;
    for (i=0;i<M;i++) {
       I[i]=i;
    }

    for (i=1; i<M; i++) {
        R = X[i];
        r = i;
		for (j=i-1; (j>=0) && (X[j]>R); j--) {
	  		X[j+1] = X[j];
            I[j+1] = I[j];
        }
		X[j+1] = R;
        I[j+1] = r;
    }
    for (i=0; i<M; i++) {
      S[I[i]]=i;
    }
    return S;
}

//data normalisation (uniform)
function  uniform(data) {
  var temp=[];
  var I = insertionSort(data);
  
  for (var i=0; i<data.length; i++) {
    temp.push(I[i]/data.length*Math.PI);
  }
  return temp;
}

//return approx. p-values from standard normal distribution
function getPval(Tval) {
	z=Tval;
	var p=1+ z*(0.04986735+ z*(0.02114101+ z*(0.00327763+ z*(0.0000380036+ z*(0.0000488906+ z*0.000005383)))));
	p=p*p; p=p*p; p=p*p;
	return (1/(p*p))/2;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//some other useful functions
////////////////////////////////////////////////////////////////////////////////////////////////

//get q quantile of the data (VaR approximation)
function getQuantile(q,data) {
var sorted = data.slice().sort(function (a, b) { return a - b; }); //sort the data
	if (q.length) {
		return quantileSorted(sorted, q);
	}
	else {
		return null;
	}
}

//sort data to get quantiles
function quantileSorted(sample, p) {
        var idx = (sample.length) * p;
        if (p < 0 || p > 1) {
            return null;
        } else if (p == 1) {
			return sample[sample.length - 1];
        } else if (p == 0) {
            return sample[0];
        } else if (idx % 1 !== 0) {
            return sample[Math.ceil(idx) - 1];
        } else if (sample.length % 2 == 0) {
            return (sample[idx - 1] + sample[idx]) / 2;
        } else {
            return sample[idx];
        }
}

//returns the absolute values of the distance to a point of multivariate vector
function abs(point,data) {
	var distance=[];
	for (var i = 0; i < data.length; i++) {
        distance.push(Math.abs(point-data[i]));
	}
    return distance;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// some simple math functions
//////////////////////////////////////////////////////////////////////////////////////////////// 

//mean
function mean(data) {
        if (data.length === 0) return null;
        return sum(data)/data.length;
}

//standard deviation
function sd(data) {
		if (data.length === 0) return null;
        var mean_value = mean(data);
        var deviations = [];
 
        for (var i=0; i<data.length; i++) {
            deviations.push(Math.pow(data[i] - mean_value, 2));
        }
        return Math.sqrt(mean(deviations));
}

//variance
function myVar(data) {
		if (data.length === 0) return null;
        var mean_value = mean(data);
        var deviations = [];
        for (var i=0; i<data.length; i++) {
            deviations.push(Math.pow(data[i] - mean_value, 2));
        }
        return mean(deviations);
}

//summation
function sum(data) {
    var value = 0;
    for (var i = 0; i < data.length; i++) {
        value += parseFloat(data[i]);
    }
    return value;
}

//returns the maximum of an array
function maximum(data) {
    var value;
    for (var i = 0; i < data.length; i++) {
        if (data[i] > value || value == undefined) value = data[i];
    }
    return value;
}
