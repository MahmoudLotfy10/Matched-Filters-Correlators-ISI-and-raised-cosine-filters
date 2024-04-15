%%%%%%%%%%%%%%%%%%% Requirement 1 %%%%%%%%%%%%%%%%%%%
ts = 1;        
Rs = 1/ts;    
p = [5 4 3 2 1]/sqrt(55);
n_bits = 10;
n_samples = 5;
data = randi ([0,1],1,n_bits);
polar_data = (2*data)-1;
data_upsampled = upsample(polar_data,n_samples);
y = conv(data_upsampled, p);

figure;
subplot(4,1,1);
plot(p);
grid on;
title("pulse shaping function");   
xlabel("Time");
ylabel("Amplitude");

subplot(4,1,2);
plot(data_upsampled);
grid on;
title("Data Upsampled");   
xlabel("Time");
ylabel("Amplitude");

subplot(4,1,3);
plot(y);
grid on;
title("Transmitted signal"); 
xlabel("Time");
ylabel("Amplitude");

%%%%%%%%%%%%%%%%%%%%  Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matched_filter = fliplr(p);
non_matched_filter = [1 1 1 1 1]/sqrt(5);

%%%%%%%%%%%%%%%%%%%%% output of filters%%%%%%%%%%%%%%%%%%%%%%
matched_filter_out = conv(y,matched_filter);
non_matched_filter_out = conv(y,non_matched_filter);
%%%%%%%%%%%%%%%%%%%%%Creation of Train Impulses%%%%%%%%%%%%%%
train_impulses = repmat([0 0 0 0 1],1,floor( size(matched_filter_out,2)/n_samples) );
subplot(4,1,4);
stem(train_impulses);
grid on;
title("train impulses"); 
xlabel("Time");
ylabel("Amplitude");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to adjust length
matched_filter_out = matched_filter_out(1:length(train_impulses));
non_matched_filter_out = non_matched_filter_out(1:length(train_impulses));

matched_filter_out_s     = matched_filter_out .* train_impulses ;
non_matched_filter_out_s = non_matched_filter_out .* train_impulses ;

correlator_out = zeros(1, length(y) - length(p) + 6);
for i = 1:5:length(y) - length(p) + 1
    out = 0;
    for j = 1:length(p)
        out = out + y(i+j-1) * p(j);
        correlator_out(i+j-1) = out; 
    end
end
corr_out_s = correlator_out .* train_impulses ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Filters');
subplot(2,1,1);
plot(matched_filter,'r');
grid on
xlabel('time');
ylabel('amplitude');
title('matched filter');
subplot(2,1,2);
plot(non_matched_filter,'b');
grid on
xlabel('time');
ylabel('amplitude');
title('non-matched filter');


figure('Name','Comparison between the outputs of the filters at the sampling instants');
subplot(2,1,1);
plot(matched_filter_out,'r');
hold on;
stem(matched_filter_out,'r');
hold off;
grid on
xlabel('time');
ylabel('amplitude');
title('matched filter output');
subplot(2,1,2);
plot(non_matched_filter_out,'b');
hold on;
stem(non_matched_filter_out,'b');
hold off;
grid on
xlabel('time');
ylabel('amplitude');
title('non-matched filter output');

figure('Name',' Comparison between the outputs of the filters after sampling');
plot(matched_filter_out_s,'r-');
hold on;
plot(non_matched_filter_out_s,'b--');
hold off;
grid on;
xlabel('time');
ylabel('amplitude');
title('Comparison between the outputs of the filters after sampling ');
legend 'matched filter' 'non matched filter';

figure('Name','output of matched and correlator');
plot(matched_filter_out,'r');
hold on;
plot(correlator_out,'b');
hold off;
grid on;
xlabel('time');
ylabel('amplitude');
title('output of matched and correlator');
legend 'matched filter' 'correlator';

figure('Name','output of matched and correlator');
subplot(2,1,1);
plot(matched_filter_out,'r');
hold on;
stem(matched_filter_out,'r')
hold off;
grid on;
xlabel('time');
ylabel('amplitude');
title('output of matched ');

subplot(2,1,2);
plot(correlator_out,'b');
hold on;
stem(correlator_out,'b')
hold off;
grid on;
xlabel('time');
ylabel('amplitude');
title('output of correlator ');

figure('Name','output of matched and correlator after sampling');
plot(matched_filter_out_s,'r-');
hold on;
plot(corr_out_s,'b--');
hold off;
grid on;
xlabel('time');
ylabel('amplitude');
title('output of matched and correlator after sampling');
legend 'matched filter' correlator;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Noise Analysis && Requirement 2 %%%%%%%%%%%%%%%%%%%

n_bits = 10^4;
n_samples = 5;
data = randi ([0,1],[1,n_bits]);
polar_data = (2*data)-1;
data_upsampled = upsample(polar_data,n_samples);

p = [5 4 3 2 1]/sqrt(55);

y = conv(data_upsampled, p);
y = y(y~=0);

unity_noise = randn(1,5*n_bits);

matched_filter = fliplr(p);
non_matched_filter_2 = [5 5 5 5 5]/sqrt(125);

snr = -2:1:5;

BER_theoretical = zeros(size(snr));
BER_1 = zeros(size(snr));
BER_2 = zeros(size(snr));

for i=snr
    index=i+3;
    no=1/db2mag(i);
    mean = 0;
    var=sqrt(no/2);
    noise = mean + var *unity_noise;
    v= y + noise;
    
    v_matched_filter    = conv(v,matched_filter);
    v_non_matched_filter= conv(v,non_matched_filter_2);
    
    v_matched_filter_out = zeros(1, size(polar_data, 2));
    v_non_matched_filter_out = zeros(1, size(polar_data, 2));
    
    for j = 1:1:size(polar_data, 2)
       if v_matched_filter(5*j) > 0
           v_matched_filter_out(j) = 1;
       else
           v_matched_filter_out(j) = -1;
       end
       
       if v_non_matched_filter(5*j) > 0
           v_non_matched_filter_out(j) = 1;
       else
           v_non_matched_filter_out(j) = -1;
       end 
   end
    temp = v_matched_filter_out ~= polar_data;
    BER_1(index) = sum(temp)/n_bits;
    
    temp = v_non_matched_filter_out ~= polar_data;
    BER_2(index) = sum(temp)/n_bits;
    
    BER_theoretical(index) = 0.5*erfc(sqrt(1/no));
    
end

figure("name" , 'BER plot');
subplot(2,2,1);
semilogy(snr,BER_1); 
title("BER matched filter");   
grid on;
xlabel("Eb/No");
ylabel("BER");

subplot(2,2,2); 
semilogy (snr,BER_2);
title("BER non matched filter");   
grid on;
xlabel("Eb/No");
ylabel("BER");

subplot(2,2,3);
semilogy(snr,BER_theoretical);
title("BER theoretical"); 
grid on;
xlabel("Eb/No");
ylabel("BER");

subplot(2,2,4);
semilogy(snr,BER_theoretical,'b');
hold on;
semilogy(snr,BER_1,'r'); 
semilogy (snr,BER_2,'g');
hold off;
title("BER Comparision"); 
legend 'BER theoretical' 'BER matched filter' 'BER non matched filter';
grid on;
xlabel("Eb/No");
ylabel("BER");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ISI && Requirement 3 %%%%%%%%%%%%%%%%%%%

n_bits = 100; 
data = randi([0, 1], [1, n_bits]);
polar_data = (2 * data) - 1;
data_upsampled = upsample(polar_data, n_samples);

R_cases = [0, 0, 1, 1];
delay_cases = [2, 8, 2, 8];

figure("Name", 'rcosine filter');
for i = 1:4 
    R = R_cases(i);
    delay = delay_cases(i);
    
    % the square root raised cosine filter
    filter = rcosine(1, n_samples, 'sqrt', R, delay);
    figure(9);
    subplot(2, 2, i);
    plot(filter);
    title("rcosine filter R: " + R + " Delay: " + delay);   
    xlabel("time samples");
    ylabel("Amplitude");
        
    transmit_data = conv(data_upsampled, filter);
    receive_data = conv(transmit_data, filter);  
    
    % Ensure transmit_data and receive_data have the same length
    min_length = min(length(transmit_data), length(receive_data));
    transmit_data = transmit_data(1:min_length);
    receive_data = receive_data(1:min_length);
    
    eye_diagram_data=[transmit_data; receive_data];
    eye_diagram=eyediagram(eye_diagram_data', n_samples * 2);
    set(eye_diagram,'Name',"Eye Diagram when R="+R_cases(i)+" & delay="+delay_cases(i)+"");
end
