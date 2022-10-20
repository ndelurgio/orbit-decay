gsurf = 980.665;
re    = 6356.77;

sav = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,...
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

% CSW
sw = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,...
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
swc = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,...
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

pt = [ ...
	% pt1[50] =
	 0.996040e+00, 0.385528e-01, 0.303445e-02,-0.105531e+00,-0.607134e-02,...
	-0.516278e-03,-0.115622e+00, 0.202240e-02, 0.990156e-02,-0.127371e+00,...
	-0.302449e-01, 0.123512e-01,-0.526277e-02,-0.845398e+01, 0.000000e+00,...
	 0.142370e-01, 0.000000e+00, 0.125818e+03, 0.805486e-02, 0.164419e-02,...
	-0.621452e-05, 0.311701e-02, 0.000000e+00, 0.386578e-02, 0.132397e+00,...
	 0.213315e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.641110e-05,...
	 0.000000e+00, 0.300150e+02, 0.533297e-02, 0.389146e-02, 0.204725e-02,...
	 0.000000e+00, 0.000000e+00,-0.192645e-01, 0.275905e+01, 0.147284e-02,...
	 0.341345e-03,-0.117388e-02,-0.354589e-03, 0.113139e+00, 0.169134e+00,...
	 0.508295e-02, 0.365016e-04, 0.426385e-02, 0.115102e-03, 0.511819e-02,...% pt2[50] =
	 0.609108e-02, 0.404995e-04, 0.153049e-02, 0.241470e-04, 0.230764e-02,...
	 0.155267e-02, 0.133722e-02,-0.182318e-02,-0.263007e+03, 0.000000e+00,...
	 0.137337e-02, 0.995774e-03, 0.000000e+00,-0.108983e+03, 0.562606e-02,...
	 0.594053e-02, 0.109358e-02, 0.000000e+00,-0.133410e-01,-0.243409e-01,...
	-0.135688e-01, 0.311370e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.283023e+04, 0.845583e-03, 0.538706e-03, 0.000000e+00, 0.247956e+03,...
	 0.292246e-02, 0.000000e+00, 0.000000e+00, 0.747703e-04, 0.887993e-03,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.116540e-01,-0.449173e-02,-0.353189e-03,-0.173933e-03,-0.153218e-03,...
	-0.565411e+00, 0.777272e-02,-0.911784e+02, 0.645187e-03, 0.000000e+00,...% pt3[50] =
	-0.837685e-03, 0.242318e-02, 0.473796e-02,-0.301801e-02,-0.423564e-02,...
	-0.248289e-02, 0.919286e-03, 0.216372e-02, 0.863968e-03, 0.189689e-02,...
	 0.415654e-02, 0.000000e+00, 0.118068e-01, 0.331190e-02, 0.000000e+00,...
	 0.120222e-02, 0.000000e+00, 0.000000e+00,-0.307246e+01, 0.000000e+00,...
	 0.000000e+00, 0.672403e-03, 0.108930e-02, 0.972278e-03, 0.468242e+01,...
	-0.315034e-03, 0.400059e-02, 0.515036e-02, 0.162989e-02, 0.108824e-02,...
	 0.995261e-03, 0.418955e+01,-0.364059e+00, 0.170182e-02, 0.000000e+00,...
	 0.000000e+00,-0.320120e+01, 0.000000e+00, 0.580206e-02, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00];
pd = [...
	% pa1[50] =
	 0.104934e+01,-0.288362e-01,-0.207095e+00,-0.103314e+00,-0.702373e-02,...
	 0.129664e-01, 0.408853e+00,-0.919895e-02,-0.188660e-01, 0.140927e+01,...
	 0.175033e+00, 0.187351e-01, 0.110979e+00,-0.742871e+01, 0.000000e+00,...
	 0.267143e+00,-0.595979e-01, 0.105038e+03,-0.840963e-01,-0.697632e-03,...
	 0.206521e-05, 0.765306e-03, 0.000000e+00, 0.000000e+00, 0.126762e+00,...
	 0.128876e+00,-0.504479e-01,-0.130735e-01,-0.224348e-01, 0.000000e+00,...
	 0.000000e+00,-0.150832e+03,-0.629928e-02, 0.000000e+00,-0.407760e-02,...
	 0.000000e+00, 0.000000e+00, 0.525725e-01,-0.311486e+02,-0.313351e-02,...
	 0.275838e-02, 0.000000e+00, 0.000000e+00, 0.111247e+00, 0.108815e+00,...
	-0.466713e-01, 0.000000e+00,-0.329329e-02, 0.000000e+00, 0.167838e-02,...% pa2[50] =	
	-0.916691e-02, 0.345044e-04,-0.971806e-02, 0.000000e+00,-0.204672e-02,...
	-0.786899e-02,-0.798285e-02, 0.536515e-02,-0.531172e+04, 0.000000e+00,...
	-0.642781e-02,-0.171690e-02, 0.000000e+00,-0.679131e+02,-0.179912e-01,...
	-0.158305e-01,-0.712313e-02, 0.000000e+00, 0.253477e-01, 0.852960e-01,...
	 0.102163e+00, 0.295009e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.684625e+04,-0.619098e-02,-0.269289e-02, 0.000000e+00,-0.520231e+03,...
	-0.633463e-02, 0.000000e+00, 0.000000e+00,-0.602428e-02,-0.407077e-02,...
	 0.542264e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.407560e-01, 0.282288e-01, 0.908088e-02, 0.000000e+00, 0.000000e+00,...
	-0.405204e+00,-0.597931e-01,-0.731823e+02,-0.206620e-02, 0.000000e+00,...% pa3[50] =
	-0.372723e-02,-0.188146e-01,-0.101794e-01, 0.804633e-02, 0.101090e-01,...
	 0.873253e-02, 0.238268e-01, 0.480444e-02, 0.171088e-02, 0.396369e-01,...
	-0.213809e-01, 0.000000e+00,-0.102588e+00,-0.591702e-02, 0.000000e+00,...
	 0.270923e-02, 0.000000e+00, 0.000000e+00,-0.175043e+03, 0.603489e+00,...
	-0.617589e+00, 0.838098e-02, 0.183871e-02,-0.705329e-03,-0.406644e+01,...
	-0.509347e-02,-0.284344e-01,-0.124160e-01, 0.133665e-01, 0.393410e-02,...
	-0.503723e-03,-0.457683e+01,-0.529542e+00,-0.425812e-02, 0.000000e+00,...
	 0.000000e+00, 0.191541e+02, 0.000000e+00, 0.323247e-02, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00;   % pb1[50] =
	 0.931113e+00,-0.138721e+00,-0.133457e+00,-0.529542e-01,-0.444983e-02,...
	 0.135264e-01, 0.598075e-01,-0.362880e-01,-0.312798e-01, 0.372068e+00,...
	 0.295974e-01, 0.120509e-01, 0.521995e-01,-0.778888e+01, 0.000000e+00,...
	 0.118634e+00,-0.204495e-01, 0.103280e+03, 0.982432e-01, 0.477694e-03,...
	 0.000000e+00, 0.274372e-02, 0.000000e+00, 0.000000e+00, 0.757809e-01,...
	 0.171403e+00,-0.105205e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00,-0.873348e+01,-0.581094e-02, 0.000000e+00,-0.814944e-02,...
	 0.000000e+00, 0.000000e+00, 0.517255e-01,-0.153028e+02,-0.348932e-02,...
	 0.961771e-03, 0.557732e-02,-0.454180e-03, 0.988213e-01, 0.940456e-01,...
	-0.318797e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.232122e-02,...% pb2[50] =
	-0.600220e-02, 0.277654e-04,-0.322019e-02, 0.000000e+00,-0.378551e-02,...
	-0.334809e-02,-0.170668e-02, 0.000000e+00, 0.636184e+04, 0.000000e+00,...
	 0.159986e-02,-0.388204e-02,-0.164825e-02,-0.747955e+02,-0.105360e-01,...
	-0.945723e-02,-0.159824e-02,-0.706730e-03,-0.168513e-01,-0.113023e+00,...
	-0.636637e-01,-0.137709e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.152368e+05,-0.586061e-02,-0.253108e-02, 0.000000e+00,-0.254837e+04,...
	-0.328988e-02, 0.000000e+00, 0.000000e+00,-0.276364e-02, 0.967923e-02,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.434255e-01, 0.114020e-01,-0.618447e-02, 0.000000e+00, 0.000000e+00,...
	-0.302568e+00,-0.327694e-01,-0.671589e+02,-0.228340e-02, 0.000000e+00,...% pb3[50] =
	 0.306230e-02,-0.465113e-02,-0.973421e-02, 0.128326e-01, 0.788553e-02,...
	 0.797197e-02,-0.120760e-01,-0.767547e-02,-0.120755e-02,-0.298523e-01,...
	-0.126560e-01, 0.000000e+00,-0.568350e-01,-0.153039e-01, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.242911e-02,-0.401347e-02,-0.219074e-02, 0.311281e+01,...
	 0.323251e-02,-0.639523e-02,-0.663069e-02,-0.304403e-03,-0.401920e-02,...
	-0.118708e-02, 0.415211e+01,-0.201896e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00;   % pc1[50] =
	 0.106903e+01, 0.377113e-03, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.898481e-01,-0.236325e+02, 0.208180e-01, 0.139638e+03,-0.119444e+00,...
	-0.845398e+01,-0.399776e-05, 0.000000e+00, 0.366210e-02,-0.178929e-02,...
	 0.190412e-01,-0.392257e-01, 0.632343e-02, 0.548144e-02, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.243022e-02,...
	 0.976619e+00, 0.568478e-03, 0.582026e-02, 0.000000e+00, 0.621998e-02,...
	 0.000000e+00, 0.000000e+00, 0.107674e-01, 0.893820e+02,-0.192414e-01,...
	-0.845398e+01, 0.000000e+00, 0.000000e+00,-0.200200e-01,-0.195833e-02,...
	-0.938391e-02, 0.131480e-01,-0.260147e-02,-0.808556e-03, 0.511651e-04,...
	 0.255717e-02, 0.000000e+00, 0.466814e-02, 0.664196e-02, 0.000000e+00,...% pc2[50] =
	 0.998594e+00, 0.190038e-03, 0.000000e+00,-0.243825e-01, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.522105e-01,...
	-0.845398e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.767271e-02, 0.564539e-02,-0.270623e-02,-0.526454e-03, 0.137075e-02,...
	 0.133060e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.949197e+00, 0.000000e+00, 0.000000e+00,-0.768008e-01, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00,-0.137993e-01,-0.140136e+01, 0.120481e+00,...
	-0.845398e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.987746e-02, 0.175330e-02,-0.688835e-03, 0.287022e-02, 0.000000e+00,...
	 0.000000e+00, 0.744513e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...% pc3[50] =
	 0.152840e+00, 0.000000e+00, 0.000000e+00, 0.116252e+01, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.649190e+00,...
	-0.845398e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.584949e-01,-0.102105e+00, 0.299153e-01,-0.486227e-01, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00;   % pd1[50] =
	 0.931402e+00, 0.137976e+00, 0.000000e+00, 0.323736e-03, 0.000000e+00,...
	-0.910906e-02, 0.707506e-01, 0.000000e+00,-0.516650e-01, 0.689755e-01,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.845398e+01, 0.000000e+00,...
	 0.281140e-01, 0.000000e+00, 0.736009e+02, 0.596604e-01, 0.000000e+00,...
	 0.000000e+00,-0.151792e-02, 0.000000e+00, 0.000000e+00, 0.132397e+00,...
	 0.213315e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.948758e+01, 0.884541e-02, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.113139e+00, 0.169134e+00,...
	 0.145192e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...% pd2[50] =
	 0.107906e-01, 0.299942e-04, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.148930e-01,...
	-0.787184e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.683420e-01,-0.441778e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.229730e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...% pd3[50] =
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00;	 % pe1[50] =
	 0.868053e+00, 0.236364e+00, 0.134306e+00, 0.103086e-01, 0.000000e+00,...
	-0.379164e-02,-0.157806e+00, 0.000000e+00,-0.587644e-01,-0.312508e+00,...
	 0.000000e+00, 0.437387e-01,-0.354091e-01,-0.223636e+02, 0.000000e+00,...
	-0.533976e-01, 0.000000e+00, 0.114091e+03, 0.517497e-01, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.132397e+00,...
	 0.213315e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.342702e+03, 0.157033e-01, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.366278e-02,...
	-0.116193e-02, 0.000000e+00, 0.000000e+00, 0.113139e+00, 0.169134e+00,...
	 0.178431e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...% pe2[50] =
	 0.162864e-01, 0.316963e-04, 0.127968e-01, 0.000000e+00, 0.000000e+00,...
	-0.704599e-02, 0.207921e-02, 0.636660e-02, 0.229940e+05, 0.000000e+00,...
	 0.127833e-01,-0.208036e-02,-0.461820e-02,-0.629391e+02,-0.120745e-01,...
	 0.136675e-01, 0.136011e-01,-0.537162e-02, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.192509e+05, 0.835522e-02, 0.419439e-02, 0.000000e+00, 0.120366e+05,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.100034e-01,-0.233267e-02,...
	 0.972374e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.265079e-01,-0.209125e-01,-0.109465e-01, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.217252e-01,-0.712385e+02,-0.189428e-02, 0.000000e+00,...% pe3[50] =
	-0.602006e-02, 0.169058e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.290646e-01,...
	 0.348971e-02, 0.000000e+00, 0.501174e-01, 0.550595e-01, 0.000000e+00,...
	-0.955897e-02, 0.000000e+00, 0.000000e+00,-0.151693e+04, 0.000000e+00,...
	 0.000000e+00, 0.129306e-01, 0.269567e-02, 0.000000e+00, 0.392243e+01,...
	-0.847690e-02, 0.116896e-01, 0.000000e+00, 0.148967e-01, 0.544521e-02,...
	 0.000000e+00, 0.564918e+01, 0.000000e+00,-0.772178e-02, 0.000000e+00,...
	 0.000000e+00,-0.734042e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00;	 % pf1[50] =
	 0.127515e+01,-0.210472e+00,-0.177924e+00, 0.218900e+00, 0.288436e-01,...
	 0.190077e-01, 0.291001e+00, 0.217437e-01,-0.105186e-01, 0.436141e+00,...
	 0.107605e+00, 0.330755e-01, 0.400581e-01,-0.958051e+01, 0.000000e+00,...
	 0.154028e-01, 0.000000e+00, 0.734194e+02, 0.496540e-01,-0.595906e-02,...
	 0.384512e-04,-0.136000e-01, 0.000000e+00, 0.000000e+00, 0.132397e+00,...
	 0.213315e+00,-0.416610e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.146276e+03,-0.198408e-01, 0.000000e+00, 0.132530e-01,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.104687e-03,...
	-0.147562e-02, 0.000000e+00, 0.000000e+00, 0.113139e+00, 0.169134e+00,...
	-0.126913e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.608370e-02,...% pf2[50] =
	-0.257587e-01, 0.319022e-04, 0.000000e+00, 0.000000e+00, 0.156644e-01,...
	 0.103640e-01, 0.105771e-02, 0.000000e+00, 0.357949e+04, 0.000000e+00,...
	-0.125672e-02, 0.152783e-02, 0.130518e-02, 0.755558e+01,-0.920341e-02,...
	-0.209142e-01,-0.134106e-01, 0.000000e+00,-0.483312e-01, 0.830900e-01,...
	 0.988009e-01,-0.141148e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.105513e+04, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.673442e-02, 0.201691e-02,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.598019e-01, 0.633298e-02,-0.112871e-02, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00,-0.128604e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...% pf3[50] =
	-0.494960e-02,-0.136415e-01,-0.115039e-01, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00,-0.586860e-02,-0.141732e-02, 0.213697e-02, 0.263845e+01,...
	-0.834186e-02,-0.187336e-01,-0.190870e-01,-0.803810e-02,-0.284279e-02,...
	 0.256722e-02, 0.171429e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00;	 % pg1[50] =
	 0.573587e+02,-0.398747e+00, 0.000000e+00,-0.529554e+00,-0.582186e-02,...
	 0.714177e-01,-0.679279e+00,-0.167715e+00,-0.642434e-01,-0.211569e+00,...
	-0.159922e+00,-0.171024e-03,-0.115885e+00, 0.651603e+01, 0.000000e+00,...
	-0.176683e+00, 0.650395e-01, 0.143504e+01, 0.928208e-01, 0.511662e-02,...
	 0.000000e+00, 0.995121e-02, 0.000000e+00, 0.000000e+00, 0.132397e+00,...
	 0.213315e+00, 0.101451e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.567667e+02, 0.238192e-02, 0.000000e+00,-0.188240e-01,...
	 0.000000e+00, 0.000000e+00, 0.476218e-01, 0.235206e+02, 0.475901e-02,...
	 0.576162e-02, 0.151815e-01,-0.192730e-01, 0.113139e+00, 0.169134e+00,...
	-0.288771e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.118418e-02,...% pg2[50] =
	-0.368927e-02, 0.314704e-04, 0.882198e-02, 0.000000e+00,-0.192562e-01,...
	-0.258674e-02,-0.219913e-01, 0.000000e+00, 0.438655e+04, 0.000000e+00,...
	 0.760126e-02, 0.259438e-02, 0.172310e-02, 0.779204e+02, 0.797786e-03,...
	-0.770510e-02, 0.190982e-02, 0.272707e-02, 0.101016e-01, 0.116537e+00,...
	-0.312236e-02, 0.139783e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	-0.130712e+04, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.320544e-02,-0.206970e-01,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.159010e-01,-0.191427e-02,-0.342829e-01, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00,-0.345379e-01, 0.894518e+02, 0.171556e-02, 0.000000e+00,...% pg3[50] =
	-0.765278e-02,-0.208987e-03,-0.157393e-01, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00,-0.860673e-02,-0.119922e-01,-0.646356e-02,-0.300107e+01,...
	-0.932511e-02,-0.150205e-01,-0.867835e-02,-0.764801e-02,-0.131495e-01,...
	-0.676720e-02,-0.182396e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00];
ps = [...
	% ph1[50] =
	 0.951363e+00,-0.467542e-01, 0.120260e+00, 0.000000e+00, 0.000000e+00,...
	 0.191357e-01, 0.000000e+00, 0.000000e+00, 0.125429e-02,-0.133240e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.845398e+01, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.252317e-02, 0.000000e+00,-0.973404e-02, 0.132397e+00,...
	 0.213315e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.718482e-03, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.787683e-02,-0.233698e-02, 0.113139e+00, 0.169134e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...% ph2[50] =
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...% ph3[50] =
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
	 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00];
pdl_0 = [...
	% pi1[50] =
	 0.933804e+00, 0.547446e+01, 0.153263e+00, 0.919303e+00, 0.164109e+02,...
	 0.427083e+01];
pdl_1 = [...
	% pi1[50] =
	 0.115897e+01, 0.471094e+00, 0.109459e+01, 0.525012e+01, 0.100000e+01,...
	 0.100000e+01, 0.103999e+01, 0.767132e+00, 0.110514e+01, 0.175636e+01,...
	 0.110845e+01, 0.233439e+01, 0.796532e+00, 0.431520e+01, 0.407300e+01,...
	 0.101885e+01, 0.239547e+00, 0.253791e-05, 0.842931e+00, 0.104192e+01,...
	 0.200202e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01];
ptm = [...
	 0.104130e+04, 0.386000e+03, 0.190000e+03, 0.166728e+02, 0.115000e+03,...
	 0.120000e+03, 0.945537e+02, 0.000000e+00];
pdm = [...
	0.245600e+08, 0.859400e+11, 0.281000e+12, 0.330000e+11, 0.133000e+10,...
	0.176100e+06, 0.100000e+08;
    0.671072e-05, 0.540000e+00, 0.000000e+00, 0.268270e+00, 0.119615e-01,...
    0.100000e+01, 0.100000e+01;
    0.100000e+03, 0.105000e+03, 0.105000e+03, 0.105000e+03, 0.105000e+03,...
    0.950000e+02, 0.105000e+03;
    0.000000e+00,-0.800000e+01, 0.280000e+02, 0.000000e+00, 0.000000e+00,...
    -0.800000e+01,-0.800000e+01;
    0.110000e+03, 0.110000e+03, 0.289500e+02, 0.110000e+03, 0.110000e+03,...
    0.110000e+03, 0.110000e+03;
	0.100000e+02, 0.100000e+02, 0.000000e+00, 0.100000e+02, 0.100000e+02,...
	0.100000e+02, 0.100000e+02;
    0.000000e+00, 0.900000e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
    0.900000e+02, 0.900000e+02;
    0.000000e+00, 0.200000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
    0.200000e+01, 0.200000e+01];
