function grad = make_derivatives_2D_complete(Size)
% Sparse matrix approximation of gradient operator on image plane.
% Use forward differences inside image, backward differences at left/bottom border

   row = Size(1);col = Size(2);
   long = row*col;
   [x,y] = meshgrid(1:col,1:row);

   linIdx = reshape(sub2ind(Size,y,x),fliplr(Size))';
   i = [reshape(linIdx(:,1:col-1).',[],1);reshape(linIdx(:,1:col-1).',[],1);
       linIdx(:,col);linIdx(:,col)];
   j = [reshape(linIdx(:,1:col-1).',[],1);reshape(linIdx(:,2:col).',[],1);
       linIdx(:,col);linIdx(:,col-1)];
   v = [ones(row*(col - 1),1)*(-1);ones(row*(col - 1),1);
       ones(row,1);ones(row,1)*(-1)];
   grad.Kx = sparse(i,j,v,long,long);
   
   i = [reshape(linIdx(1:row-1,:).',[],1);reshape(linIdx(1:row-1,:).',[],1);
       reshape(linIdx(row,:).',[],1);reshape(linIdx(row,:).',[],1)];
   j = [reshape(linIdx(1:row-1,:).',[],1);reshape(linIdx(2:row,:).',[],1);
       reshape(linIdx(row,:).',[],1);reshape(linIdx(row - 1,:).',[],1)];
   v = [ones((row-1)*col,1)*(-1);ones((row-1)*col,1);
       ones(col,1);ones(col,1)*(-1)];
   grad.Ky = sparse(i,j,v,long,long);
end