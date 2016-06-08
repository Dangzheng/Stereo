function result = huber_function(data, epsilon)
   idx1 = abs(data) >= epsilon; 
   idx2 = abs(data) < epsilon;
   result = zeros(size(data));
   result(idx1) = abs(data(idx1)) - 0.5*epsilon;
   result(idx2) = (abs(data(idx2)).^2)/(2*epsilon);
   
end