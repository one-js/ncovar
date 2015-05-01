////////////////////////////////////////////////////////////////////////////////////////////////
//
// ncovar-2.0.js
//
// JavaScript function to calculate the test statistic and p-values of the NCoVaR test from Diks and Wolski (2014)
// 
// main function: ncovar(dataA,dataB,VaR,safe,bandwidth,forwStart,forward,backward,scenario,variance,kernel,norm) { }
// dataA: time series for variable A
// dataB: time series for variable B (has to be the same length)
// VaR: Value-at-Risk level (for instance 0.05)
// safe: safe quantile level (for instance 0.5)
// bandwidth: bandwidth value
// forwStart: beginning of the forward-looking variable (0 for contemporary)
// forw: end of the forward-looking variable (0 for contemporary)
// backward: number of lags (0 for contemporary)
// scenario: "sc" for scenario 1 (tail to tail) and "sc2" for scenario 2 (median to tail)
// variance: variance estimator ("standard" or "HAC" for heteroskedasticity and autocorrelation consistent estimator)
// kernel: kernel function, default is the "rectangular" kernel (other kernels: "gaussian", "triangular", "epanechnikov", "quartic", "triweight", "tricube", "cosine", "logistic"
// norm: data normalization ("standard" for standard formal and "uniform" for uniform transformation)
//
////////////////////////////////////////////////////////////////////////////////////////////////

//main function
function ncovar(dataA,dataB,VaR,safe,bandwidth,forwStart,forward,backward,scenario,variance,kernel,norm) {
//data normalisation

if(norm=="standard") {
	var A=normalise(dataA.slice(0,dataA.length-forward));
	var B=normalise(dataB.slice(0,dataB.length-forward));
	var As=normalise(dataA.slice());
	var Bs=normalise(dataB.slice());
}
else {
	var A=uniform(dataA.slice(0,dataA.length-forward));
	var B=uniform(dataB.slice(0,dataB.length-forward));
	var As=uniform(dataA.slice());
	var Bs=uniform(dataB.slice());
}

//get quantiles depending on scenario considered
var qrA=getQuantile(VaR,A); //qxir
var qsA=getQuantile(safe,A); //qxis
var qrB=getQuantile(VaR,normalise(dataB.slice(forward,dataB.length))); //qz
if(scenario=="sc1") {
	var qsB=getQuantile(VaR,B); //qxj
}
else if(scenario=="sc2") {
	var qsB=getQuantile(safe,B); //qxj
}

//calculate T-values
return(tVal(As,Bs,forwStart,forward,backward,qrB,qsB,qrA,qsA,bandwidth,variance,kernel));
}

//calculate test statistic
function tVal(dataA,dataB,forwStart,forw,back,qz,qxj,qxir,qxis,eps,variance,kernel) {
var n=dataA.length;
var forL=forw-forwStart+1;
var K1=[];
for(var i = 0; i < n; i++) {
    K1.push(0.0);
}

for(var k=back; k<=(n-forw); k++) {
	var CIzr=0, CIzs=0, CIr=0, CIs=0;
	var Izr=0, Izs=0, Ir=0, Is=0;
	
	//calculate kernels (I) for index k
	switch(kernel) {
		case "gaussian":
			Ir=gaussianKernel(qxj,dataB.slice(k-back,k+1),eps)*gaussianKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*gaussianKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=gaussianKernel(qxj,dataB.slice(k-back,k+1),eps)*gaussianKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*gaussianKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
		break;
		case "triangular":
			Ir=trKernel(qxj,dataB.slice(k-back,k+1),eps)*trKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*trKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=trKernel(qxj,dataB.slice(k-back,k+1),eps)*trKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*trKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
		break;
		case "epanechnikov":
			Ir=epKernel(qxj,dataB.slice(k-back,k+1),eps)*epKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*epKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=epKernel(qxj,dataB.slice(k-back,k+1),eps)*epKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*epKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
		break;
		case "quartic":
			Ir=qKernel(qxj,dataB.slice(k-back,k+1),eps)*qKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*qKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=qKernel(qxj,dataB.slice(k-back,k+1),eps)*qKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*qKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
		break;
		case "triweight":
			Ir=triweightKernel(qxj,dataB.slice(k-back,k+1),eps)*triweightKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*triweightKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=triweightKernel(qxj,dataB.slice(k-back,k+1),eps)*triweightKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*triweightKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
		break;
		case "tricube":
			Ir=tricubeKernel(qxj,dataB.slice(k-back,k+1),eps)*tricubeKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*tricubeKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=tricubeKernel(qxj,dataB.slice(k-back,k+1),eps)*tricubeKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*tricubeKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
		break;
		case "cosine":
			Ir=cosKernel(qxj,dataB.slice(k-back,k+1),eps)*cosKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*cosKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=cosKernel(qxj,dataB.slice(k-back,k+1),eps)*cosKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*cosKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
		break;
		case "logistic":
			Ir=logKernel(qxj,dataB.slice(k-back,k+1),eps)*logKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*logKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=logKernel(qxj,dataB.slice(k-back,k+1),eps)*logKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*logKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
		break;
		default: 
			Ir=normKernel(qxj,dataB.slice(k-back,k+1),eps)*normKernel(qxir,dataA.slice(k-back,k+1),eps);
			Izr=Ir*normKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
			Is=normKernel(qxj,dataB.slice(k-back,k+1),eps)*normKernel(qxis,dataA.slice(k-back,k+1),eps);
			Izs=Is*normKernel(qz,dataB.slice(k+forwStart,k+forw+1),eps);
	}
	
	var m=back;
	while (m<k) {
		//calculate correlation kernels (CI) for index m
		switch(kernel) {
			case "gaussian":
				CIr=CIr+gaussianKernel(qxj,dataB.slice(m-back,m+1),eps)*gaussianKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+gaussianKernel(qxj,dataB.slice(m-back,m+1),eps)*gaussianKernel(qxir,dataA.slice(m-back,m+1),eps)*gaussianKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+gaussianKernel(qxj,dataB.slice(m-back,m+1),eps)*gaussianKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+gaussianKernel(qxj,dataB.slice(m-back,m+1),eps)*gaussianKernel(qxis,dataA.slice(m-back,m+1),eps)*gaussianKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
			break;
			case "triangular":
				CIr=CIr+trKernel(qxj,dataB.slice(m-back,m+1),eps)*trKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+trKernel(qxj,dataB.slice(m-back,m+1),eps)*trKernel(qxir,dataA.slice(m-back,m+1),eps)*trKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+trKernel(qxj,dataB.slice(m-back,m+1),eps)*trKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+trKernel(qxj,dataB.slice(m-back,m+1),eps)*trKernel(qxis,dataA.slice(m-back,m+1),eps)*trKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
			break;
			case "epanechnikov":
				CIr=CIr+epKernel(qxj,dataB.slice(m-back,m+1),eps)*epKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+epKernel(qxj,dataB.slice(m-back,m+1),eps)*epKernel(qxir,dataA.slice(m-back,m+1),eps)*epKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+epKernel(qxj,dataB.slice(m-back,m+1),eps)*epKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+epKernel(qxj,dataB.slice(m-back,m+1),eps)*epKernel(qxis,dataA.slice(m-back,m+1),eps)*epKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
			break;
			case "quartic":
				CIr=CIr+qKernel(qxj,dataB.slice(m-back,m+1),eps)*qKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+qKernel(qxj,dataB.slice(m-back,m+1),eps)*qKernel(qxir,dataA.slice(m-back,m+1),eps)*qKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+qKernel(qxj,dataB.slice(m-back,m+1),eps)*qKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+qKernel(qxj,dataB.slice(m-back,m+1),eps)*qKernel(qxis,dataA.slice(m-back,m+1),eps)*qKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
			break;
			case "triweight":
				CIr=CIr+triweightKernel(qxj,dataB.slice(m-back,m+1),eps)*triweightKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+triweightKernel(qxj,dataB.slice(m-back,m+1),eps)*triweightKernel(qxir,dataA.slice(m-back,m+1),eps)*triweightKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+triweightKernel(qxj,dataB.slice(m-back,m+1),eps)*triweightKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+triweightKernel(qxj,dataB.slice(m-back,m+1),eps)*triweightKernel(qxis,dataA.slice(m-back,m+1),eps)*triweightKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
			break;
			case "tricube":
				CIr=CIr+tricubeKernel(qxj,dataB.slice(m-back,m+1),eps)*tricubeKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+tricubeKernel(qxj,dataB.slice(m-back,m+1),eps)*tricubeKernel(qxir,dataA.slice(m-back,m+1),eps)*tricubeKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+tricubeKernel(qxj,dataB.slice(m-back,m+1),eps)*tricubeKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+tricubeKernel(qxj,dataB.slice(m-back,m+1),eps)*tricubeKernel(qxis,dataA.slice(m-back,m+1),eps)*tricubeKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
			break;
			case "cosine":
				CIr=CIr+cosKernel(qxj,dataB.slice(m-back,m+1),eps)*cosKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+cosKernel(qxj,dataB.slice(m-back,m+1),eps)*cosKernel(qxir,dataA.slice(m-back,m+1),eps)*cosKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+cosKernel(qxj,dataB.slice(m-back,m+1),eps)*cosKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+cosKernel(qxj,dataB.slice(m-back,m+1),eps)*cosKernel(qxis,dataA.slice(m-back,m+1),eps)*cosKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
			break;
			case "logistic":
				CIr=CIr+logKernel(qxj,dataB.slice(m-back,m+1),eps)*logKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+logKernel(qxj,dataB.slice(m-back,m+1),eps)*logKernel(qxir,dataA.slice(m-back,m+1),eps)*logKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+logKernel(qxj,dataB.slice(m-back,m+1),eps)*logKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+logKernel(qxj,dataB.slice(m-back,m+1),eps)*logKernel(qxis,dataA.slice(m-back,m+1),eps)*logKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
			break;
			default:
				CIr=CIr+normKernel(qxj,dataB.slice(m-back,m+1),eps)*normKernel(qxir,dataA.slice(m-back,m+1),eps);
				CIzr=CIzr+normKernel(qxj,dataB.slice(m-back,m+1),eps)*normKernel(qxir,dataA.slice(m-back,m+1),eps)*normKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
				CIs=CIs+normKernel(qxj,dataB.slice(m-back,m+1),eps)*normKernel(qxis,dataA.slice(m-back,m+1),eps);
				CIzs=CIzs+normKernel(qxj,dataB.slice(m-back,m+1),eps)*normKernel(qxis,dataA.slice(m-back,m+1),eps)*normKernel(qz,dataB.slice(m+forwStart,m+forw+1),eps);
		
		}
		m++;	
	}
	//Correlation Integral
	K1[k]=(Izr*CIs-Ir*CIzs+CIzr*Is-CIr*Izs);
	K1[k]=Math.pow(2*eps,(-forL-2*(back+1)-2*(back+1)))/Math.pow(n-forL,2)*K1[k];
}
var T=sum(K1.slice(0,n-forL));
//variance of the test statistics
if(variance=="HAC") {
	var tVar=HACvar(K1.slice(0,n-forL));
}
else {
	var tVar=myVarK1(K1.slice(0,n-forL));
}
var T_VAL=T*Math.sqrt(n-forL)/Math.sqrt(4*tVar);
return T_VAL;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//kernel functions
////////////////////////////////////////////////////////////////////////////////////////////////

//Gaussian
function gaussianKernel(point,vector,eps) {
var value=1;
	for(var i=0;i<vector.length;i++) {
		var p=(point-vector[i])/eps;
		value=value*(1/((Math.sqrt(2*Math.PI)*eps))*Math.exp(-1/2*Math.pow(p,2)));
	}
return value;
}

//rectangular
function normKernel(point,vector,eps) {
var value=1;
	for(var i=0; i<vector.length; i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(Math.abs(p)<1) {
			control=1;
		}
		value=Number(value*control);
	}
return value;
}

//triweight
function trKernel(point,vector,eps) {
var value=1;
	for(var i=0;i<vector.length;i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(Math.abs(p)<1) {
			control=1;
		}
		value=value*(1-Math.abs(p))*control;
	}
return value;
}

//Epanechnikov
function epKernel(point,vector,eps) {
var value=1;
	for(var i=0;i<vector.length;i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(Math.abs(p)<1) {
			control=1;
		}
		value=value*3/4*(1-Math.pow(p,2))*control;
	}
return value;
}

//quartic (biweight)
function qKernel(point,vector,eps) {
var value=1;
	for(var i=0;i<vector.length;i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(Math.abs(p)<1) {
			control=1;
		}
		value=value*15/16*Math.pow((1-Math.pow(p,2)),2)*control;
	}
return value;
}

//triweight
function triweightKernel(point,vector,eps) {
var value=1;
	for(var i=0;i<vector.length;i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(Math.abs(p)<1) {
			control=1;
		}
		value=value*35/32*Math.pow((1-Math.pow(p,2)),3)*control;
	}
return value;
}

//tricube
function tricubeKernel(point,vector,eps) {
var value=1;
	for(var i=0;i<vector.length;i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(Math.abs(p)<1) {
			control=1;
		}
		value=value*70/81*Math.pow((1-Math.pow(Math.abs(p),3)),3)*control;
	}
return value;
}

//cosine
function cosKernel(point,vector,eps) {
var value=1;
	for(var i=0;i<vector.length;i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(Math.abs(p)<1) {
			control=1;
		}
		value=value*Math.PI/4*Math.cos(Math.PI*p)*control;
	}
return value;
}

//logistic
function logKernel(point,vector,eps) {
var value=1;
	for(var i=0;i<vector.length;i++) {
		var p=(Math.abs(point-vector[i]))/eps;
		var control=0;
		if(p<1) {
			control=1;
		}
		value=value*(1/(Math.exp(p)+2+Math.exp(-p)));
	}
return value;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// data normalization functions
////////////////////////////////////////////////////////////////////////////////////////////////

//data normalization (standard normal)
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

//data normalization (uniform)
function  uniform(data) {
  var temp=[];
  var I = insertionSort(data);
  
  for (var i=0; i<data.length; i++) {
    temp.push(I[i]/data.length*Math.PI);
  }
  return temp;
}

//HAC estimator of variance
function HACvar(K1) {
	if (K1.length === 0) return null;
	n=K1.length;
	var out=0;
	var B=Math.floor(Math.pow(n,(1/4)));
	var R=[];
	var w=[];
	
	//empty matrices
	for(var i=0; i<B; i++) {
		R[i] = 0;
		w[i] = 0;
	}
	
	w[0]=1;
	for(var b=0; b<B; b++) {
	//weighting function
	if(b>0) {
		w[b]=2-(2*b/B);
	}
	//sample covariance over the first index
	R[b]=0;
	for (i=b;i!=n;i++) {
        R[b] += K1[i]*K1[i-b]*(n*n); //scaled as K1 is multipled by 1/n
    }
    R[b] /= (n-b);
	}
	  	
	for(var i=0; i<B; i++) {
		out+=w[i]*R[i];
	}
    return out;
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
var tempData = data.slice();
var sorted = tempData.sort(function(a, b){return a-b;}); //sort the data
var r = 0.0;
	if (q.length) {
		r = quantileSorted(sorted, q);
	}
return r;
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
        //} else if (idx % 1 !== 0) {
            //return sample[Math.ceil(idx) - 1];
        //} else if (sample.length % 2 == 0) {
            //return (sample[idx - 1] + sample[idx]) / 2;
        } else {
        	var h = (p*((sample.length)-1)+1)-1; //negative one as the counter starts at 0
			var hfloor = Math.floor(h);
			var quant = sample[hfloor]+(h-hfloor)*(sample[hfloor+1]-sample[hfloor]);
	    	return quant;
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
//K1 variance
function myVarK1(data) {
		if (data.length === 0) return null;
		var n=data.length;
        var mean_value = mean(data);
        var deviations = [];
        for (var i=0; i<n; i++) {
            deviations.push(Math.pow(data[i] - mean_value, 2)*n*n);
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
