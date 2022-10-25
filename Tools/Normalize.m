function OutImg = Normalize(InImg)
	% all element divide the max
	A=tril(InImg);
	for i=1:size(A,2)
		A(i:end,i)=A(i:end,i)/A(i,i);
	end
	C=tril(A,-1);
	OutImg=C'+A;
end
