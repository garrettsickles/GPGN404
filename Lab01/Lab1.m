% NOTE: In all cases the domain was an array of integers from -20 to 20 (inclusive)
function Lab1(domain)
exp8(domain);
end
% Problem 1: Impulse ?[n]
function impulse(domain)
val = zeros(size(domain));
for i=1:size(domain)
    if domain(i) == 0
        val(i) = 1;
    end
end
stem(domain,val)
end


% Problem 2:  Unit Step u[n]
function unitStep(domain)
val = zeros(size(domain));
for i=1:size(domain)
    if domain(i) >= 0
        val(i) = 1;
    else
        val(i) < 0;
    end
end
stem(domain,val)
end

% Problem 3: Exponential Decay x[n] = (a^n)*u[n] with a = 0.9
function exponentialDecay(domain)
val = zeros(size(domain));
for i=1:size(domain)
    if domain(i) >= 0
        val(i) = 1;
    else
        val(i) < 0;
    end
    val(i) = val(i)*((0.9)^domain(i));
end
stem(domain,val)
end

% Problem 4: Sine for ?_0 = ?/20.0 and ? = 0
function sin4(domain)
val = zeros(size(domain));
for i=1:size(domain)
    val(i) = sin((3.14159/20.0)*domain(i)+0.0);
end
stem(domain,val)
end

% Problem 5: Cosine for ?_0 = ?/20.0 and ? = 0
function cos5(domain)
val = zeros(size(domain));
for i=1:size(domain)
    val(i) = cos((3.14159/20.0)*domain(i)+0.0);
end
stem(domain,val)
end

% Problem 6: Complex Exponential for ?_0 = ?/20.0 and ? = 0
function exp6(domain)
re = zeros(size(domain));
im = zeros(size(domain));
for j=1:size(domain)
    val = exp(((3.14159/20.0)*domain(j)+0.0)*1i);
    re(j) = real(val);
    im(j) = imag(val);
end
stem(domain, re, 'color', 'b'); hold on;
stem(domain, im, 'color', 'r');
end

% Problem 7: Complex Exponential for ?_0 = ?/20.0 and ? = ??/2.0
function exp7(domain)
re = zeros(size(domain));
im = zeros(size(domain));
for j=1:size(domain)
    val = exp(((3.14159/20.0)*domain(j)-(3.14159/2.0))*1i);
    re(j) = real(val);
    im(j) = imag(val);
end
stem(domain, re, 'color', 'b'); hold on;
stem(domain, im, 'color', 'r');
end

% Problem 8: Complex Exponential for ? 0 = ?/10.0 and ? = 0.0
function exp8(domain)
re = zeros(size(domain));
im = zeros(size(domain));
for j=1:size(domain)
    val = exp(((3.14159/10.0)*domain(j)+0.0)*1i);
    re(j) = real(val);
    im(j) = imag(val);
end
stem(domain, re, 'color', 'b'); hold on;
stem(domain, im, 'color', 'r');
end

% Problem 9: Sine for ?_0 = ??2/20.0 and ? = 0.0
function sin9(domain)
val = zeros(size(domain));
for j=1:size(domain)
    val(j) = sin((3.14159/20.0)*sqrt(2)*domain(j)+0.0);
end
stem(domain,val)
end

% Problem 10: Triangular discrete-time wave x[n] = mod(n + 10, 21) ? 10
function tdtw(domain)
val = zeros(size(domain));
for j=1:size(domain)
    val(j) = mod(domain(j)+10,21)-10.0;
end
stem(domain,val)
end